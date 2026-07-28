using Test
using JuMP
using QMBCertify

@testset "Phase 0 baseline and model-specific Pauli symmetry" begin
    heisenberg = heisenberg_symmetry()
    ising = ising_chain_symmetry()

    @test reduce!(UInt16[1], L=4, symmetry=heisenberg)[2] == 0
    @test reduce!(UInt16[1], L=4, symmetry=ising) == (UInt16[1], 1)
    @test reduce!(UInt16[3], L=4, symmetry=ising)[2] == 0
    @test reduce!(UInt16[3, 6], L=4, symmetry=ising)[2] != 0
    @test reduce!(UInt16[4], L=4, symmetry=ising)[1] == UInt16[1]

    reduced = symmetry_reduce_support([[1], [4], [3], [3, 6]], 4, symmetry=ising)
    @test reduced == [UInt16[1], UInt16[3, 6]]
    @test fieldnames(IsingMomentResult) == (:energy, :energy_density, :status,
        :basis_size, :moment_count, :raw_entry_count)
    @test hasmethod(GSB, Tuple{Vector{Vector{Int}}, Vector{Float64}, Int, Int})
    @test validate_relaxation_label(:basis_policy, :heuristic_basis) == :heuristic_basis
    @test validate_relaxation_label(:reduction, :equivalent_reduction) == :equivalent_reduction
    @test validate_relaxation_label(:strengthening, :valid_strengthening) == :valid_strengthening
    @test validate_relaxation_label(:result_scope, :numerical_diagnostic) == :numerical_diagnostic
    @test_throws ArgumentError validate_relaxation_label(:result_scope, :strict_certificate)
    @test_throws ArgumentError validate_relaxation_label(:unknown, :heuristic_basis)
end

@testset "Phase 0 canonical Pauli algebra and UInt16 boundary" begin
    @test pauli_product(Int[]) == (UInt16[], 1)
    @test pauli_product([1, 1]) == (UInt16[], 1)
    @test pauli_product([1, 2]) == (UInt16[3], im)
    @test pauli_product([2, 1]) == (UInt16[3], -im)
    @test pauli_product([4, 1]) == (UInt16[1, 4], 1)
    @test pauli_product([1, 2]; realify=true) == (UInt16[3], 1)
    @test pauli_product([2, 1]; realify=true) == (UInt16[3], -1)
    @test pauli_product([typemax(UInt16)]) == (UInt16[typemax(UInt16)], 1)
    @test_throws ArgumentError pauli_product([Int(typemax(UInt16)) + 1])
    @test_throws ArgumentError pauli_product([0])
end

@testset "Phase 1 explicit polynomial, specification, and symmetry declarations" begin
    polynomial = PauliPolynomial([
        UInt16[4] => 2.0,
        UInt16[1] => 1.0,
        UInt16[1] => -0.25,
        UInt16[2, 1] => im,
        UInt16[4] => -2.0,
    ])
    @test length(polynomial.terms) == 2
    @test polynomial.terms[1].word == UInt16[1]
    @test polynomial.terms[1].coefficient == 0.75
    @test polynomial.terms[2].word == UInt16[3]
    @test polynomial.terms[2].coefficient == 1
    @test_throws ArgumentError PauliPolynomial([UInt16[1] => im])
    @test_throws ArgumentError PauliPolynomial([UInt16[1] => Inf])

    sector = BasisSector(:explicit, [UInt16[], UInt16[1], UInt16[2]])
    @test sector.words == [UInt16[], UInt16[1], UInt16[2]]
    @test_throws ArgumentError BasisSector(:duplicate, [UInt16[1], UInt16[1]])
    metadata = TranslationOrbitMetadata(2, [1, 1, 2])
    @test BasisSector(:orbits, sector.words; translation_orbits=metadata).translation_orbits.period == 2

    hamiltonian = PauliPolynomial([UInt16[1] => 1.0])
    bad = SymmetryDeclaration(:bad_flip, :moment_zero, [1]; axis_sign=(-1, 1, 1))
    @test_throws ArgumentError RelaxationSpecification(hamiltonian, [sector]; symmetries=[bad])
    @test_throws ArgumentError RelaxationSpecification(hamiltonian, [sector];
        certificate_scope=:rigorously_postvalidated)
    @test_throws ArgumentError RelaxationSpecification(hamiltonian, [sector];
        certificate_scope=:unknown)
    identity_on_one_site = SymmetryDeclaration(:identity, :moment_equality, [1])
    @test_throws ArgumentError RelaxationSpecification(hamiltonian, [sector];
        symmetries=[identity_on_one_site],
        observables=Dict(:site_two => PauliPolynomial([UInt16[4] => 1.0])))
    @test_throws ArgumentError RelaxationSpecification(hamiltonian, [sector];
        observables=Dict(:nonhermitian => PauliPolynomial([UInt16[1] => im]; hermitian=false)))
end

@testset "Phase 1 deterministic B†B, strengthening, RDM, and closure" begin
    hamiltonian = PauliPolynomial([UInt16[1] => 0.7, UInt16[3, 6] => -1.2])
    basis = [BasisSector(:main, [UInt16[], UInt16[1], UInt16[2], UInt16[3]])]
    linear_tests = [
        PauliPolynomial([UInt16[2] => 1.0]),
        PauliPolynomial([UInt16[2] => 2.0]),
        hamiltonian,
    ]
    spec = RelaxationSpecification(hamiltonian, basis;
        linear_tests=linear_tests,
        psd_state_basis=[UInt16[], UInt16[2]],
        rdm_regions=[RDMRegion(:sites_1_2, [1, 2])],
        observables=Dict(:magnetization => PauliPolynomial([UInt16[1] => 1.0])))

    first_compile = compile_relaxation(spec)
    second_compile = compile_relaxation(spec)
    diagnostics = first_compile.diagnostics
    @test first_compile.moments == second_compile.moments
    @test diagnostics.fingerprint == second_compile.diagnostics.fingerprint
    changed_observable_spec = RelaxationSpecification(hamiltonian, basis;
        linear_tests=linear_tests,
        psd_state_basis=[UInt16[], UInt16[2]],
        rdm_regions=[RDMRegion(:sites_1_2, [1, 2])],
        observables=Dict(:magnetization => PauliPolynomial([UInt16[1] => 2.0])))
    @test compile_relaxation(changed_observable_spec).diagnostics.fingerprint != diagnostics.fingerprint
    @test length(diagnostics.fingerprint) == 16
    @test diagnostics.raw_basis_size == 4
    @test diagnostics.raw_bdagb_entries == 16
    main_block = first_compile.psd_blocks[1]
    @test main_block.entries[2, 3].words == [UInt16[3]]
    @test main_block.entries[2, 3].coefficients == ComplexF64[im]
    @test main_block.entries[3, 2].coefficients == ComplexF64[-im]
    @test first_compile.moment_index[UInt16[3]] == second_compile.moment_index[UInt16[3]]
    @test diagnostics.linear_rows_raw == 3
    @test diagnostics.linear_rows_duplicate == 1
    @test diagnostics.linear_rows_zero == 1
    @test diagnostics.linear_rows_independent == 1
    commutator_spec = RelaxationSpecification(
        PauliPolynomial([UInt16[1] => 0.7]),
        [BasisSector(:commutator_basis, [UInt16[], UInt16[1], UInt16[2], UInt16[3]])];
        linear_tests=[PauliPolynomial([UInt16[2] => 1.0])])
    commutator_row = only(compile_relaxation(commutator_spec).linear_rows)
    @test commutator_row.words == [UInt16[3]]
    @test commutator_row.coefficients ≈ ComplexF64[1.4im]
    @test diagnostics.psd_block_count == 3
    @test diagnostics.psd_block_dimensions == [4, 2, 4]
    @test diagnostics.max_psd_block_dimension == 4
    @test diagnostics.real_embedded_dimensions == [8, 4, 8]
    @test diagnostics.state_optimality_moment_increment >= 0
    pso_spec = RelaxationSpecification(
        PauliPolynomial([UInt16[1] => 1.0]),
        [BasisSector(:pso_moment, [UInt16[], UInt16[1], UInt16[2], UInt16[3]])];
        psd_state_basis=[UInt16[], UInt16[2]])
    pso_block = only(filter(block -> block.role == :state_optimality,
                            compile_relaxation(pso_spec).psd_blocks))
    @test isempty(pso_block.entries[1, 1].words)
    @test pso_block.entries[1, 2].words == [UInt16[3]]
    @test pso_block.entries[1, 2].coefficients == ComplexF64[im]
    @test pso_block.entries[2, 1].coefficients == ComplexF64[-im]
    @test pso_block.entries[2, 2].words == [UInt16[1]]
    @test pso_block.entries[2, 2].coefficients == ComplexF64[-2]
    @test eltype(pso_block.entries[2, 2].coefficients) == ComplexF64
    @test diagnostics.rdm_moment_increment > 0
    @test isempty(diagnostics.missing_support)
    @test diagnostics.scope == :numerical_relaxation
    @test haskey(first_compile.observables, :magnetization)
    @test any(contains("B†B"), first_compile.psd_blocks[1].provenance)
    @test any(contains("RDM"), diagnostics.constraint_provenance)

    for block in first_compile.psd_blocks, row in axes(block.entries, 1), column in axes(block.entries, 2)
        left = Dict(zip(block.entries[row, column].words, block.entries[row, column].coefficients))
        right = Dict(word => conj(coefficient) for (word, coefficient) in
                     zip(block.entries[column, row].words, block.entries[column, row].coefficients))
        @test left == right
    end
    rdm = only(filter(block -> block.role == :rdm, first_compile.psd_blocks))
    @test UInt16[] in rdm.entries[1, 1].words
    identity_trace = sum(begin
        location = findfirst(==(UInt16[]), rdm.entries[i, i].words)
        location === nothing ? 0.0 + 0.0im : rdm.entries[i, i].coefficients[location]
    end for i in axes(rdm.entries, 1))
    @test identity_trace == 1

    blocked_rdm_spec = RelaxationSpecification(hamiltonian, basis;
        rdm_regions=[RDMRegion(:blocked_site, [1]; blocks=[[1], [2]])])
    blocked_rdm = filter(block -> block.role == :rdm,
                         compile_relaxation(blocked_rdm_spec).psd_blocks)
    @test getfield.(blocked_rdm, :name) == [:blocked_site_block_1, :blocked_site_block_2]
    @test size.(getfield.(blocked_rdm, :entries), 1) == [1, 1]
    bad_partition = RelaxationSpecification(hamiltonian, basis;
        rdm_regions=[RDMRegion(:bad_partition, [1]; blocks=[[1]])])
    @test_throws ArgumentError compile_relaxation(bad_partition)

    built = build_jump_model(first_compile)
    @test built.compiled === first_compile
    @test string(JuMP.termination_status(built.model)) == "OPTIMIZE_NOT_CALLED"
    observable_built = build_jump_model(first_compile; objective=:magnetization)
    @test observable_built.compiled === first_compile
    @test string(JuMP.termination_status(observable_built.model)) == "OPTIMIZE_NOT_CALLED"
    @test_throws ArgumentError build_jump_model(first_compile; objective=:missing)

    incomplete = RelaxationSpecification(hamiltonian, basis;
        declared_moment_support=[UInt16[]])
    closure_error = try
        compile_relaxation(incomplete)
        nothing
    catch error
        error
    end
    @test closure_error isa ArgumentError
    @test contains(string(closure_error), "moment closure failed")
    @test contains(string(closure_error), "required by")
    @test contains(string(closure_error), "Hamiltonian objective")
    @test contains(string(closure_error), "B†B main")
end

@testset "Phase 1 explicit symmetry purposes and Fourier blocks" begin
    block_hamiltonian = PauliPolynomial([UInt16[3] => 1.0])
    block_generator = SymmetryDeclaration(:xy_parity, :basis_block, [1];
                                          axis_sign=(-1, -1, 1))
    block_spec = RelaxationSpecification(block_hamiltonian,
        [BasisSector(:characters, [UInt16[], UInt16[1], UInt16[3]])];
        symmetries=[block_generator])
    blocked = compile_relaxation(block_spec)
    @test blocked.diagnostics.psd_block_dimensions == [1, 2]
    @test Set(getfield.(blocked.psd_blocks, :name)) ==
          Set([:characters_character_1, Symbol("characters_character_-1")])

    translation = SymmetryDeclaration(:translation_by_one, :fourier_orbit, [2, 1])
    fourier_hamiltonian = PauliPolynomial([UInt16[1] => 1.0, UInt16[4] => 1.0])
    fourier_sector = BasisSector(:translated_x, [UInt16[1], UInt16[4]];
        translation_orbits=TranslationOrbitMetadata(2, [1, 1]))
    fourier = compile_relaxation(RelaxationSpecification(fourier_hamiltonian,
        [fourier_sector]; symmetries=[translation]))
    @test fourier.diagnostics.psd_block_dimensions == [1, 1]
    @test getfield.(fourier.psd_blocks, :name) ==
          [:translated_x_momentum_0, :translated_x_momentum_1]
    @test all(contains("Fourier"), fourier.diagnostics.constraint_provenance)

    incomplete_orbit = BasisSector(:incomplete, [UInt16[1]];
        translation_orbits=TranslationOrbitMetadata(2, [1]))
    @test_throws ArgumentError compile_relaxation(RelaxationSpecification(
        fourier_hamiltonian, [incomplete_orbit]; symmetries=[translation]))
    @test_throws ArgumentError compile_relaxation(RelaxationSpecification(
        fourier_hamiltonian, [BasisSector(:no_metadata, [UInt16[1], UInt16[4]])];
        symmetries=[translation]))

    bad_block = SymmetryDeclaration(:site_swap, :basis_block, [2, 1])
    @test_throws ArgumentError compile_relaxation(RelaxationSpecification(
        fourier_hamiltonian, [BasisSector(:nondiagonal, [UInt16[1], UInt16[4]])];
        symmetries=[bad_block]))
end

@testset "Phase 1 Ising compile/build/solve split" begin
    compiled = compile_ising_relaxation(1.3, 0.7, 4; degree=1)
    @test compiled isa CompiledRelaxation
    @test compiled.diagnostics.raw_basis_size == 13
    @test compiled.diagnostics.raw_bdagb_entries == 169
    @test compiled.diagnostics.scalar_moment_count < compiled.diagnostics.raw_bdagb_entries
    @test compiled.diagnostics.psd_block_dimensions == [13]
    @test compiled.specification.certificate_scope == :numerical_relaxation
    @test string(JuMP.termination_status(build_jump_model(compiled).model)) == "OPTIMIZE_NOT_CALLED"
end

@testset "Phase 1 Heisenberg regression and Table 2: 31 is max PSD block dimension, not constraint count" begin
    L = 4
    even_basis = QMBCertify.get_basis(L, 0, 2)
    odd_basis = QMBCertify.get_basis(L, 1, 2)
    @test length(even_basis) == 3L
    @test length(odd_basis) == 3L
    @test even_basis[1:L] == [UInt16[1, 4], UInt16[4, 7], UInt16[7, 10], UInt16[1, 10]]
    @test reduce!(copy(even_basis[1]), L=L)[1] == UInt16[1, 4]

    benchmark = heisenberg_table2_benchmark(100, 4, 1)
    @test benchmark.sparse_basis_size == 12_001
    @test length(benchmark.basis.words) == 12_001
    @test benchmark.max_psd_block_dimension == 31
    @test length(benchmark.symmetry_block_upper_bounds) == 100
    @test maximum(benchmark.symmetry_block_upper_bounds) == benchmark.max_psd_block_dimension
    @test all(==(31), benchmark.symmetry_block_upper_bounds)
    @test benchmark.original_dimension == 8_127_090_301
    @test benchmark.equality_reduced_dimension == 322_029_976
    @test benchmark.max_psd_block_dimension != benchmark.sparse_basis_size
    @test_throws ArgumentError heisenberg_table2_benchmark(100, 4, 2)
end

function mosek_license_available()
    license_paths = String[]
    haskey(ENV, "MOSEKLM_LICENSE_FILE") && append!(license_paths, split(ENV["MOSEKLM_LICENSE_FILE"], ':'))
    push!(license_paths, joinpath(homedir(), "mosek", "mosek.lic"))
    return any(isfile, license_paths)
end

@testset "Phase 2.1 explicit dimerized J1-J2 Hamiltonian" begin
    L = 6
    mg = dimerized_j1j2_hamiltonian(1.0, 0.5, 0.0, L)
    mg_terms = Dict(term.word => real(term.coefficient) for term in mg.terms)
    @test length(mg.terms) == 6L
    @test mg_terms[UInt16[1, 4]] == 0.25
    @test mg_terms[UInt16[1, 7]] == 0.125
    @test all(axis -> haskey(mg_terms, UInt16[axis, 3 + axis]), 1:3)

    decoupled = dimerized_j1j2_hamiltonian(1.0, 0.0, 1.0, L)
    decoupled_terms = Dict(term.word => real(term.coefficient) for term in decoupled.terms)
    pauli_label(site, axis) = UInt16(3 * (site - 1) + axis)
    @test length(decoupled.terms) == 3L ÷ 2
    @test all(decoupled_terms[sort(UInt16[pauli_label(site, axis),
                                                   pauli_label(mod1(site + 1, L), axis)])] == 0.5
              for site in 2:2:L, axis in 1:3)
    @test all(!haskey(decoupled_terms, sort(UInt16[pauli_label(site, axis),
                                                        pauli_label(mod1(site + 1, L), axis)]))
              for site in 1:2:(L - 1), axis in 1:3)
    @test_throws ArgumentError dimerized_j1j2_hamiltonian(1.0, 0.5, 0.0, 5)
    @test_throws ArgumentError dimerized_j1j2_hamiltonian(Inf, 0.5, 0.0, L)
end

@testset "Phase 2.2 fixed-budget uniform and dimer-adapted bases" begin
    L = 6
    budget = 1 + 12L
    uniform = dimerized_chain_basis(:uniform_local, L; budget=budget)
    adapted = dimerized_chain_basis(:operator_adapted, L; budget=budget)
    @test length(uniform.words) == length(adapted.words) == budget
    @test uniform.words[1] == adapted.words[1] == UInt16[]
    @test all(length(word) <= 2 for word in uniform.words)

    strong_x = UInt16[4, 7]
    weak_x = UInt16[1, 4]
    j2_x = UInt16[1, 7]
    adjacent_dimers_x = UInt16[4, 7, 10, 13]
    @test all(word -> word in adapted.words, (strong_x, weak_x, j2_x, adjacent_dimers_x))
    @test adjacent_dimers_x ∉ uniform.words
    @test_throws ArgumentError dimerized_chain_basis(:unknown, L; budget=budget)
    @test_throws ArgumentError dimerized_chain_basis(:operator_adapted, L; budget=10_000)
end

@testset "Phase 2.3 parameter-matched symmetry and strengthening" begin
    L = 6
    uniform_symmetries = dimerized_chain_symmetries(0.0, L)
    dimerized_symmetries = dimerized_chain_symmetries(0.4, L)
    @test uniform_symmetries[1].name == :translation_by_1
    @test uniform_symmetries[1].site_map == [2, 3, 4, 5, 6, 1]
    @test dimerized_symmetries[1].name == :translation_by_2
    @test dimerized_symmetries[1].site_map == [3, 4, 5, 6, 1, 2]

    baseline = dimerized_chain_specification(1.0, 0.5, 0.0, L;
        policy=:uniform_local, budget=25, strengthening=:baseline)
    linear = dimerized_chain_specification(1.0, 0.5, 0.0, L;
        policy=:uniform_local, budget=25, strengthening=:linear)
    psd = dimerized_chain_specification(1.0, 0.5, 0.0, L;
        policy=:uniform_local, budget=25, strengthening=:psd)
    @test isempty(baseline.linear_tests) && isempty(baseline.psd_state_basis)
    @test length(linear.linear_tests) == 1 && isempty(linear.psd_state_basis)
    @test length(psd.linear_tests) == 1 && length(psd.psd_state_basis) == 4
    @test baseline.normalization == L
    @test compile_relaxation(psd).diagnostics.psd_block_count == 2

    wrong_translation = SymmetryDeclaration(:wrong_translation, :moment_equality,
        [2, 3, 4, 5, 6, 1])
    dimerized_hamiltonian = dimerized_j1j2_hamiltonian(1.0, 0.0, 0.4, L)
    @test_throws ArgumentError RelaxationSpecification(dimerized_hamiltonian,
        [dimerized_chain_basis(:uniform_local, L; budget=25)];
        symmetries=[wrong_translation])
    @test_throws ArgumentError dimerized_chain_specification(1.0, 0.0, 0.4, L;
        policy=:uniform_local, budget=25, strengthening=:unknown)
end

@testset "license-aware solver regressions" begin
    if !mosek_license_available()
        @test_skip "Mosek license unavailable: solver-dependent Ising and GSB regressions skipped"
    else
        L = 4
        h = 1.7
        field_only = ising_ground_state_bound(0.0, h, L; degree=1, QUIET=true)
        @test string(field_only.status) == "OPTIMAL"
        @test field_only.energy ≈ -L * abs(h) atol=1e-7
        @test field_only.energy_density ≈ -abs(h) atol=1e-7
        @test field_only.moment_count < field_only.raw_entry_count

        J = 1.3
        coupling_only = ising_ground_state_bound(J, 0.0, L; degree=1, QUIET=true)
        @test string(coupling_only.status) == "OPTIMAL"
        @test coupling_only.energy ≈ -L * abs(J) atol=1e-7
        @test coupling_only.energy_density ≈ -abs(J) atol=1e-7
        @test coupling_only.moment_count < coupling_only.raw_entry_count

        value, data = GSB([[1, 4]], [0.75], 4, 2;
            lso=false, pso=0, rdm=false, QUIET=true)
        @test isfinite(value)
        @test propertynames(data) == (:correlation1, :correlation2, :correlation3,
            :basis, :sbasis, :tsupp, :GramMat, :sGramMat, :multiplier, :moment)
        @test length(data.basis) == 2
        @test UInt16[] in data.tsupp
    end
end
