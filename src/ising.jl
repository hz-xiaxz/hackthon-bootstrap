"""Result of a symmetry-reduced transverse-field Ising moment relaxation."""
struct IsingMomentResult
    energy::Float64
    energy_density::Float64
    status
    basis_size::Int
    moment_count::Int
    raw_entry_count::Int
end

function _pauli_basis(L::Int, degree::Int)
    L > 1 || throw(ArgumentError("L must be at least 2"))
    0 <= degree <= L || throw(ArgumentError("degree must lie between 0 and L"))

    basis = Vector{UInt16}[UInt16[]]
    function add_words!(sites, next_site, sites_left)
        if sites_left == 0
            for labels in Iterators.product(ntuple(_ -> 1:3, length(sites))...)
                push!(basis, UInt16[3 * (sites[k] - 1) + labels[k] for k in eachindex(sites)])
            end
            return
        end
        for site in next_site:(L - sites_left + 1)
            push!(sites, site)
            add_words!(sites, site + 1, sites_left - 1)
            pop!(sites)
        end
    end
    for support_size in 1:degree
        add_words!(Int[], 1, support_size)
    end
    return basis
end

function _ising_moment_data(basis, L, symmetry)
    n = length(basis)
    representatives = Vector{UInt16}[UInt16[]]
    representative_index = Dict{Vector{UInt16},Int}(UInt16[] => 1)
    entry_index = zeros(Int, n, n)
    entry_coefficient = zeros(ComplexF64, n, n)

    for row in 1:n, column in row:n
        word = UInt16[reverse(basis[row]); basis[column]]
        representative, coefficient = reduce!(word, L=L, symmetry=symmetry)
        coefficient == 0 && continue
        index = get!(representative_index, representative) do
            push!(representatives, representative)
            length(representatives)
        end
        entry_index[row, column] = index
        entry_coefficient[row, column] = coefficient
        entry_index[column, row] = index
        entry_coefficient[column, row] = conj(coefficient)
    end
    return representatives, representative_index, entry_index, entry_coefficient
end

"""Build the explicit transverse-field Ising relaxation specification without solving it."""
function ising_relaxation_specification(J::Real, h::Real, L::Int; degree::Int=1,
        symmetry::PauliSymmetryModel=ising_chain_symmetry(),
        certificate_scope=:numerical_relaxation)
    basis_words = _pauli_basis(L, degree)
    sector = BasisSector(:ising_moment, basis_words)

    terms = Pair{Vector{UInt16},Float64}[]
    for site in 1:L
        next_site = mod1(site + 1, L)
        push!(terms, UInt16[3 * (site - 1) + 3, 3 * (next_site - 1) + 3] => -Float64(J))
        push!(terms, UInt16[3 * (site - 1) + 1] => -Float64(h))
    end
    hamiltonian = PauliPolynomial(terms)

    declarations = SymmetryDeclaration[]
    push!(declarations, SymmetryDeclaration(:global_ising_flip, :moment_zero, collect(1:L);
                                             axis_sign=(1, -1, -1)))
    if symmetry.translation
        push!(declarations, SymmetryDeclaration(:translation, :moment_equality,
                                                 [collect(2:L); 1]))
    end
    if symmetry.reflection
        push!(declarations, SymmetryDeclaration(:reflection, :moment_equality,
                                                 [1; collect(L:-1:2)]))
    end
    return RelaxationSpecification(hamiltonian, [sector]; symmetries=declarations,
                                   certificate_scope=certificate_scope)
end

"""Compile the explicit Ising relaxation without creating a solver model."""
function compile_ising_relaxation(J::Real, h::Real, L::Int; degree::Int=1,
                                  symmetry::PauliSymmetryModel=ising_chain_symmetry())
    return compile_relaxation(ising_relaxation_specification(J, h, L;
        degree=degree, symmetry=symmetry))
end

"""
    ising_ground_state_bound(J, h, L; degree=1, symmetry=ising_chain_symmetry(),
                             QUIET=true, mosek_setting=mosek_para())

Compute a lower bound on the ground-state energy of the periodic chain
`H = -J∑ᵢ ZᵢZᵢ₊₁ - h∑ᵢ Xᵢ` using the explicit Gate A compiler and a
symmetry-reduced Pauli moment matrix. The public arguments and result shape are
preserved; compilation, JuMP model construction, and optimization are separate.
"""
function ising_ground_state_bound(J::Real, h::Real, L::Int; degree::Int=1,
                                   symmetry::PauliSymmetryModel=ising_chain_symmetry(),
                                   QUIET::Bool=true,
                                   mosek_setting::mosek_para=mosek_para())
    compiled = compile_ising_relaxation(J, h, L; degree=degree, symmetry=symmetry)
    attributes = Pair{String,Any}[
        "MSK_DPAR_INTPNT_CO_TOL_PFEAS" => min(mosek_setting.tol_pfeas, 1e-9),
        "MSK_DPAR_INTPNT_CO_TOL_DFEAS" => min(mosek_setting.tol_dfeas, 1e-9),
        "MSK_DPAR_INTPNT_CO_TOL_REL_GAP" => min(mosek_setting.tol_relgap, 1e-9),
        "MSK_IPAR_NUM_THREADS" => mosek_setting.num_threads,
    ]
    built = build_jump_model(compiled; optimizer=Mosek.Optimizer,
                             optimizer_attributes=attributes)
    set_optimizer_attribute(built.model, MOI.Silent(), QUIET)
    solved = solve_relaxation(built)
    solved.status == MOI.OPTIMAL || error("Ising moment relaxation terminated with status $(solved.status)")
    energy = solved.objective_value
    n = compiled.diagnostics.raw_basis_size
    return IsingMomentResult(
        energy,
        energy / L,
        solved.status,
        n,
        compiled.diagnostics.scalar_moment_count,
        div(n * (n + 1), 2),
    )
end
