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

"""
    ising_ground_state_bound(J, h, L; degree=1, symmetry=ising_chain_symmetry(),
                             QUIET=true, mosek_setting=mosek_para())

Compute a lower bound on the ground-state energy of the periodic chain
`H = -J∑ᵢ ZᵢZᵢ₊₁ - h∑ᵢ Xᵢ` using a symmetry-reduced Pauli moment matrix.
`degree` is the maximum number of sites in a word in the moment basis.

Translation, reflection, and global Ising spin flip identify scalar moments,
while positivity is imposed on the unreduced basis. Since the Hamiltonian is
real, moments are restricted to their conjugation-invariant (real) sector;
the complex Hermitian moment matrix is represented by a real PSD embedding.
"""
function ising_ground_state_bound(J::Real, h::Real, L::Int; degree::Int=1,
                                   symmetry::PauliSymmetryModel=ising_chain_symmetry(),
                                   QUIET::Bool=true,
                                   mosek_setting::mosek_para=mosek_para())
    basis = _pauli_basis(L, degree)
    representatives, representative_index, entry_index, entry_coefficient =
        _ising_moment_data(basis, L, symmetry)

    model = Model(optimizer_with_attributes(
        Mosek.Optimizer,
        "MSK_DPAR_INTPNT_CO_TOL_PFEAS" => mosek_setting.tol_pfeas,
        "MSK_DPAR_INTPNT_CO_TOL_DFEAS" => mosek_setting.tol_dfeas,
        "MSK_DPAR_INTPNT_CO_TOL_REL_GAP" => mosek_setting.tol_relgap,
        "MSK_IPAR_NUM_THREADS" => mosek_setting.num_threads,
    ))
    set_optimizer_attribute(model, MOI.Silent(), QUIET)

    @variable(model, moments[1:length(representatives)])
    @constraint(model, moments[1] == 1)

    n = length(basis)
    real_part = Matrix{AffExpr}(undef, n, n)
    imag_part = Matrix{AffExpr}(undef, n, n)
    for row in 1:n, column in 1:n
        index = entry_index[row, column]
        if index == 0
            real_part[row, column] = AffExpr(0.0)
            imag_part[row, column] = AffExpr(0.0)
        else
            coefficient = entry_coefficient[row, column]
            real_part[row, column] = coefficient.re * moments[index]
            imag_part[row, column] = coefficient.im * moments[index]
        end
    end
    real_embedding = [real_part -imag_part; imag_part real_part]
    @constraint(model, Symmetric(real_embedding) in PSDCone())

    function moment_for(word)
        representative, coefficient = reduce!(UInt16.(word), L=L, symmetry=symmetry)
        coefficient == 0 && return AffExpr(0.0)
        index = get(representative_index, representative, 0)
        index == 0 && throw(ArgumentError("degree=$degree does not generate Hamiltonian moment $word"))
        isreal(coefficient) || throw(ArgumentError("Hamiltonian reduced to a non-real moment"))
        return real(coefficient) * moments[index]
    end

    x_moment = moment_for(UInt16[1])
    zz_moment = moment_for(UInt16[3, 6])
    @objective(model, Min, L * (-float(J) * zz_moment - float(h) * x_moment))
    optimize!(model)

    status = termination_status(model)
    status == MOI.OPTIMAL || error("Ising moment relaxation terminated with status $status")
    energy = objective_value(model)
    return IsingMomentResult(
        energy,
        energy / L,
        status,
        n,
        length(representatives),
        div(n * (n + 1), 2),
    )
end
