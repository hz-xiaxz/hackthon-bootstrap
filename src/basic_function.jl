mutable struct mosek_para
    tol_pfeas::Float64
    tol_dfeas::Float64
    tol_relgap::Float64
    num_threads::Int64
end

mosek_para() = mosek_para(1e-8, 1e-8, 1e-8, 0)

function get_basis(L, label, d; lattice="chain", extra=0, three_type=[1;1])
    basis = Vector{UInt16}[]
    if lattice == "square"
        tb2 = [[1;0], [0;1], [1;1], [1;-1], [2;0], [0;2], [2;1], [1;2], [2;2]]
        lb2 = 9
        if L >= 6
            lb2 = 12
            push!(tb2, [2;-1], [1;-2], [2;-2], [0;3], [1;3], [2;3], [3;3], [3;2], [3;1], [3;0]) # 7
        end
        if L >= 8
            push!(tb2, [3;-1], [3;-2], [3;-3], [2;-3], [1;-3], [0;4], [1;4], [2;4], [3;4], [4;4], [4;3], [4;2], [4;1], [4;0]) # 21
        end
        if L >= 10
            push!(tb2, [4;-1], [4;-2], [4;-3], [4;-4], [3;-4], [2;-4], [1;-4], [0;5], [1;5], [2;5], [3;5], [4;5], [5;5], [5;4], [5;3], [5;2], [5;1], [5;0]) # 39
        end
        if three_type == [1;1]
            tb3 = [[0;1;1;1], [0;1;-1;1], [1;0;1;1], [-1;0;-1;1], [1;0;2;0], [0;1;0;2]]
        else
            tb3 = [[0;1;1;1], [0;1;-1;1], [1;0;1;1], [-1;0;-1;1], [1;1;2;2], [-1;1;-2;2]]
        end
    end
    if label > 0
        a1 = [[rot(label)[1];rot(label)[2]], [rot(label)[2];rot(label)[1]]]
        a2 = [[label;label;rot(label)[1];rot(label)[2]], [rot(label)[2];rot(label)[1];label;label], [label;rot(label)[1];label;rot(label)[2]], [rot(label)[2];label;rot(label)[1];label], [label;rot(label)[1];rot(label)[2];label], [label;rot(label)[2];rot(label)[1];label], 
        [rot(label)[1];label;label;rot(label)[2]], [rot(label)[2];label;label;rot(label)[1]], [rot(label)[1];label;rot(label)[2];label], [label;rot(label)[2];label;rot(label)[1]], [rot(label)[1];rot(label)[2];label;label], [label;label;rot(label)[2];rot(label)[1]],  
        [rot(label)[1];rot(label)[1];rot(label)[1];rot(label)[2]], [rot(label)[2];rot(label)[1];rot(label)[1];rot(label)[1]], [rot(label)[1];rot(label)[1];rot(label)[2];rot(label)[1]], [rot(label)[1];rot(label)[2];rot(label)[1];rot(label)[1]], 
        [rot(label)[2];rot(label)[2];rot(label)[2];rot(label)[1]], [rot(label)[1];rot(label)[2];rot(label)[2];rot(label)[2]], [rot(label)[2];rot(label)[2];rot(label)[1];rot(label)[2]], [rot(label)[2];rot(label)[1];rot(label)[2];rot(label)[2]]]
        if lattice == "chain"
            for i = 1:L
                push!(basis, [3*(i-1)+label])
            end
            if d > 2
                for i = 1:L
                    push!(basis, sort([3*(i-1)+label;smod(3*(i-1+three_type[1])+label, 3L);smod(3*(i-1+sum(three_type))+label, 3L)]))
                end
                for l = 1:2, k = 1:3
                    ind = rot(label)[l]*ones(Int, 3)
                    ind[k] = label
                    for i = 1:L
                        push!(basis, sort([3*(i-1)+ind[1];smod(3*(i-1+three_type[1])+ind[2], 3L);smod(3*(i-1+sum(three_type))+ind[3], 3L)]))
                    end
                end
            end
            if d > 1
                for s = 0:extra, k in a1, i = 1:L
                    push!(basis, sort([3*(i-1)+k[1];smod(3*(i+s)+k[2], 3L)]))
                end
            end
            if d > 3
                for k in a2, i = 1:L
                    push!(basis, sort([3*(i-1)+k[1];smod(3*i+k[2], 3L);smod(3*(i+1)+k[3], 3L);smod(3*(i+2)+k[4], 3L)]))
                end
            end
        else
            for i = 1:L, j = 1:L
                push!(basis, [3*(slabel(i, j, L=L)-1)+label])
            end
            if d > 2
                for s in tb3, i = 1:L, j = 1:L
                    push!(basis, sort([3*(slabel(i, j, L=L)-1)+label;3*(slabel(i+s[1], j+s[2], L=L)-1)+label;3*(slabel(i+s[3], j+s[4], L=L)-1)+label]))
                end
                for k = 1:3, l = 1:2
                    ind = rot(label)[l]*ones(Int, 3)
                    ind[k] = label
                    for s in tb3, i = 1:L, j = 1:L
                        push!(basis, sort([3*(slabel(i, j, L=L)-1)+ind[1];3*(slabel(i+s[1], j+s[2], L=L)-1)+ind[2];3*(slabel(i+s[3], j+s[4], L=L)-1)+ind[3]]))
                    end
                end
            end
            if d > 1
                for k in a1, s in tb2[1:lb2+extra], i = 1:L, j = 1:L
                    push!(basis, sort([3*(slabel(i, j, L=L)-1)+k[1];3*(slabel(i+s[1], j+s[2], L=L)-1)+k[2]]))
                end
            end
            if d > 3
                for k in a2, i = 1:L, j = 1:L
                    push!(basis, sort([3*(slabel(i, j, L=L)-1)+k[1];3*(slabel(i+1, j, L=L)-1)+k[2];3*(slabel(i, j+1, L=L)-1)+k[3];3*(slabel(i+1, j+1, L=L)-1)+k[4]]))
                end
            end
        end  
    else
        a1 = [[1;2;3], [3;2;1], [1;3;2], [2;3;1], [2;1;3], [3;1;2]]
        a2 = [[1;1;1;1], [2;2;2;2], [3;3;3;3], [1;2;2;1], [2;1;1;2], [1;3;3;1], [3;1;1;3], [3;2;2;3], [2;3;3;2], [1;1;2;2], [2;2;1;1], [1;2;1;2], [2;1;2;1], [1;1;3;3], [3;3;1;1], [1;3;1;3], [3;1;3;1], [3;3;2;2], [2;2;3;3], [3;2;3;2], [2;3;2;3]]
        if lattice == "chain"
            for s = 0:extra, k = 1:3, i = 1:L
                push!(basis, sort([3*(i-1)+k;smod(3*(i+s)+k, 3L)]))
            end
            if d > 3
                for k in a2, i = 1:L
                    push!(basis, sort([3*(i-1)+k[1];smod(3*i+k[2], 3L);smod(3*(i+1)+k[3], 3L);smod(3*(i+2)+k[4], 3L)]))
                end
            end
            if d > 2
                for k in a1, i = 1:L
                    push!(basis, sort([3*(i-1)+k[1];smod(3*(i-1+three_type[1])+k[2], 3L);smod(3*(i-1+sum(three_type))+k[3], 3L)]))
                end
            end
        else
            for s in tb2[1:lb2+extra], k = 1:3, i = 1:L, j = 1:L
                push!(basis, sort([3*(slabel(i, j, L=L)-1)+k;3*(slabel(i+s[1], j+s[2], L=L)-1)+k]))
            end
            if d > 3
                for k in a2, i = 1:L, j = 1:L
                    push!(basis, sort([3*(slabel(i, j, L=L)-1)+k[1];3*(slabel(i+1, j, L=L)-1)+k[2];3*(slabel(i, j+1, L=L)-1)+k[3];3*(slabel(i+1, j+1, L=L)-1)+k[4]]))
                end
            end
            if d > 2
                for s in tb3, k in a1, i = 1:L, j = 1:L
                    push!(basis, sort([3*(slabel(i, j, L=L)-1)+k[1];3*(slabel(i+s[1], j+s[2], L=L)-1)+k[2];3*(slabel(i+s[3], j+s[4], L=L)-1)+k[3]]))
                end
            end
        end
    end
    return basis
end

# binary search in a sorted sequence
function bfind(A, a)
    low = 1
    high = length(A)
    while low <= high
        mid = ceil(Int, (low+high)/2)
        if A[mid] == a
           return mid
        elseif A[mid] < a
           low = mid + 1
        else
           high = mid - 1
        end
    end
    return nothing
end

# reduction to the normal form
function reduce1!(a::Vector{UInt16})
    la = length(a)
    flag = 1
    while flag == 1
        ind = findfirst(x->ceil(Int, a[x]/3) > ceil(Int, a[x+1]/3), 1:la-1)
        if ind !== nothing
            a[ind],a[ind+1] = a[ind+1],a[ind]
            flag = 1
        else
            flag = 0
        end
    end
    return a
end

# reduction to the normal form
function reduce2!(a::Vector{UInt16}; realify=false)
    la = length(a)
    flag = 1
    coef = 1
    while flag == 1
        ind = findfirst(x -> a[x] != a[x+1] && ceil(Int, a[x]/3) == ceil(Int, a[x+1]/3), 1:la-1)
        if ind !== nothing
            s = mod.(a[ind:ind+1], 3)
            if s == [1, 2]
                a[ind] += UInt16(2)
                coef *= im
            elseif s == [0, 2]
                a[ind] -= UInt16(2)
                coef *= -im
            elseif s == [2, 1] || s == [1, 0]
                a[ind] += UInt16(1)
                coef *= -im
            else
                a[ind] -= UInt16(1)
                coef *= im  
            end
            deleteat!(a, ind+1)
            la -= 1
            flag = 1
        else
            flag = 0
        end
    end
    if realify == true && !isreal(coef)
        coef = imag(coef)
    end
    return a,coef
end

# reduction to the normal form
function reduce3!(a::Vector{UInt16})
    i = 1
    while i < length(a)
        if a[i] == a[i+1]
            deleteat!(a, i)
            deleteat!(a, i)
        else
            i += 1
        end
    end
    return a
end

abstract type AbstractSymmetryModel end

"""
Symmetries used to identify Pauli words and eliminate moments.

`axis_permutations` acts on Pauli labels `(X,Y,Z)`. Each `sign_generator`
is a tuple whose true entries change sign under one independent Z₂ symmetry.
Spatial translations/reflections are handled separately from internal symmetry.
"""
struct PauliSymmetryModel <: AbstractSymmetryModel
    axis_permutations::Vector{NTuple{3,UInt8}}
    sign_generators::Vector{NTuple{3,Bool}}
    translation::Bool
    reflection::Bool
end

const _IDENTITY_AXIS_PERMUTATION = (UInt8(1), UInt8(2), UInt8(3))

"""Symmetry assumptions historically used by `GSB` for Heisenberg models."""
function heisenberg_symmetry(; translation=true, reflection=true)
    axis_permutations = NTuple{3,UInt8}[
        (1, 2, 3), (1, 3, 2), (2, 1, 3),
        (2, 3, 1), (3, 1, 2), (3, 2, 1),
    ]
    sign_generators = NTuple{3,Bool}[(true, true, false), (false, true, true)]
    return PauliSymmetryModel(axis_permutations, sign_generators, translation, reflection)
end

"""
Symmetry model for the periodic transverse-field Ising chain
`H = -J∑ ZᵢZᵢ₊₁ - h∑ Xᵢ`.

The global `∏Xᵢ` symmetry changes the signs of Y and Z. Consequently every
word containing an odd total number of Y/Z factors has zero invariant moment.
"""
function ising_chain_symmetry(; translation=true, reflection=true)
    return PauliSymmetryModel(
        NTuple{3,UInt8}[_IDENTITY_AXIS_PERMUTATION],
        NTuple{3,Bool}[(false, true, true)],
        translation,
        reflection,
    )
end

# Identify moments forced to zero by model-specific sign symmetries.
function isz(a::Vector{UInt16}, symmetry::PauliSymmetryModel=heisenberg_symmetry())
    counts = ntuple(label -> count(==(label), smod.(a, 3)), 3)
    return any(symmetry.sign_generators) do generator
        isodd(sum(counts[label] for label in 1:3 if generator[label]))
    end
end

function perm(a, symmetry::PauliSymmetryModel=heisenberg_symmetry())
    labels = smod.(a, 3)
    sites = 3 .* (ceil.(Int, a ./ 3) .- 1)
    return [UInt16.(sites .+ [axis_permutation[label] for label in labels])
            for axis_permutation in symmetry.axis_permutations]
end

# Reduction with respect to spatial and internal symmetries.
function reduce4(a::Vector{UInt16}, L; lattice="chain", symmetry::PauliSymmetryModel=heisenberg_symmetry())
    isempty(a) && return a
    L > 0 || throw(ArgumentError("L must be positive when reducing a nonempty word by spatial symmetry"))

    spatial_orbit = Vector{UInt16}[]
    if lattice == "chain"
        if symmetry.translation
            for i in eachindex(a)
                translated = [a[i:end]; a[1:i-1] .+ 3L] .- 3 * (ceil(UInt16, a[i] / 3) - 1)
                push!(spatial_orbit, translated)
                if symmetry.reflection
                    reversed = reverse(translated)
                    reflected = 3 .* (ceil(UInt16, translated[end] / 3) .- ceil.(UInt16, reversed / 3)) + smod.(reversed, 3)
                    push!(spatial_orbit, reflected)
                end
            end
        else
            push!(spatial_orbit, copy(a))
            if symmetry.reflection
                sites = ceil.(UInt16, a / 3)
                reflected = 3 .* (sites[end] .- reverse(sites)) + smod.(reverse(a), 3)
                push!(spatial_orbit, reflected)
            end
        end
    else
        factor = [[1;1], [-1;1], [1;-1], [-1;-1]]
        loc = location.(ceil.(UInt16, a / 3))
        anchors = symmetry.translation ? eachindex(a) : (1:1)
        transforms = symmetry.reflection ? (1:4) : (1:1)
        for i in anchors, k in transforms
            temp = zeros(UInt16, length(a))
            for j in eachindex(a)
                p = slabel(factor[k][1] * (loc[j][1] - loc[i][1]) + 1,
                           factor[k][2] * (loc[j][2] - loc[i][2]) + 1, L=L)
                temp[j] = 3p + a[j] - 3 * ceil(UInt16, a[j] / 3)
            end
            push!(spatial_orbit, sort(temp))
            if symmetry.reflection
                for j in eachindex(a)
                    p = slabel(factor[k][1] * (loc[j][2] - loc[i][2]) + 1,
                               factor[k][2] * (loc[j][1] - loc[i][1]) + 1, L=L)
                    temp[j] = 3p + a[j] - 3 * ceil(UInt16, a[j] / 3)
                end
                push!(spatial_orbit, sort(temp))
            end
        end
    end

    orbit = Vector{UInt16}[]
    for word in spatial_orbit
        append!(orbit, perm(word, symmetry))
    end
    return minimum(orbit)
end

# Implement Pauli algebra reduction followed by model-specific symmetry reduction.
function reduce!(a::Vector{UInt16}; L=0, lattice="chain", realify=false,
                 symmetry::PauliSymmetryModel=heisenberg_symmetry())
    reduce1!(a)
    reduce3!(a)
    a,coef = reduce2!(a, realify=realify)
    reduce3!(a)
    if isz(a, symmetry)
        coef = 0
    else
        a = reduce4(a, L, lattice=lattice, symmetry=symmetry)
    end
    return a,coef
end

"""Return the distinct nonzero symmetry representatives of Pauli words."""
function symmetry_reduce_support(words, L; lattice="chain", symmetry=heisenberg_symmetry())
    representatives = Vector{UInt16}[]
    for word in words
        representative, coefficient = reduce!(UInt16.(word), L=L, lattice=lattice, symmetry=symmetry)
        coefficient == 0 || push!(representatives, representative)
    end
    sort!(representatives)
    unique!(representatives)
    return representatives
end

function slabel(i, j; L=0)
    i = mod(i, L)==0 ? L : mod(i, L)
    j = mod(j, L)==0 ? L : mod(j, L)
    r = max(i,j)
    return r == i ? (r-1)^2+j : r^2+1-i
end

function location(p)
    r = ceil(Int, sqrt(p))
    if p-(r-1)^2 <= r
        return r, p-(r-1)^2
    else
        return r^2+1-p, r
    end
end

function rot(label)
    if label == 1
        return 2,3
    elseif label == 2
        return 3,1
    else
        return 1,2
    end
end

function smod(i, s)
    r = mod(i, s)
    return r == 0 ? s : r
end

# compute the eigenvalues of a (symmetry, real) circulant matrix
function eigen_circmat(supp, coe, L; symmetry=false, real_matrix=false)
    ne = real_matrix == true ? Int(L/2)+1 : L
    seig = [Vector{UInt16}[] for i = 1:ne]
    if symmetry == false
        ceig = [ComplexF64[] for i = 1:ne]
        for i = 1:ne
            for j = 1:L, (s,c) in enumerate(coe[j])
                if c != 0
                    push!(seig[i], supp[j][s])
                    push!(ceig[i], c*(cos(2*pi*(i-1)*(j-1)/L) + sin(2*pi*(i-1)*(j-1)/L)*im))
                end
            end
            if !isempty(ceig[i])
                seig[i], ceig[i] = resort(seig[i], ceig[i])
            end
        end
    else
        ceig = [Float64[] for i = 1:ne]
        for i = 1:ne
            for (s,c) in enumerate(coe[1])
                if c != 0
                    push!(seig[i], supp[1][s])
                    push!(ceig[i], c)
                end
            end
            for j = 2:Int(L/2), (s,c) in enumerate(coe[j])
                if c != 0
                    push!(seig[i], supp[j][s])
                    push!(ceig[i], 2*c*cos(2*pi*(i-1)*(j-1)/L))
                end
            end
            for (s,c) in enumerate(coe[end])
                if c != 0
                    push!(seig[i], supp[end][s])
                    push!(ceig[i], c*(-1)^(i-1))
                end
            end
            if !isempty(ceig[i])
                seig[i], ceig[i] = resort(seig[i], ceig[i])
            end
        end
    end
    return seig, ceig
end

# -----------------------------------------------------------------------------
# Explicit deterministic Pauli relaxation compiler (Gate A)
# -----------------------------------------------------------------------------

const PAULI_AXES = (UInt8(1), UInt8(2), UInt8(3))
const BASIS_POLICY_KINDS = (:heuristic_basis,)
const REDUCTION_KINDS = (:equivalent_reduction,)
const STRENGTHENING_KINDS = (:valid_strengthening,)
const CERTIFICATE_SCOPES = (:numerical_relaxation, :solver_bound, :numerical_diagnostic,
                            :rigorously_postvalidated)
const RELAXATION_TERMINOLOGY = (basis_policy=BASIS_POLICY_KINDS,
    reduction=REDUCTION_KINDS, strengthening=STRENGTHENING_KINDS,
    result_scope=CERTIFICATE_SCOPES)
const SYMMETRY_PURPOSES = (:moment_zero, :moment_equality, :basis_block,
                           :fourier_orbit, :redundant_block_equivalence)

"""Validate one of the four non-interchangeable relaxation terminology classes."""
function validate_relaxation_label(kind::Symbol, value)
    hasproperty(RELAXATION_TERMINOLOGY, kind) || throw(ArgumentError("unknown terminology class $kind"))
    label = Symbol(value)
    label in getproperty(RELAXATION_TERMINOLOGY, kind) ||
        throw(ArgumentError("invalid $kind label $label"))
    return label
end

"""Canonical Pauli-polynomial term. Words use `3(site-1)+axis`, with `1=X,2=Y,3=Z`."""
struct PauliTerm
    word::Vector{UInt16}
    coefficient::ComplexF64
end

_word_key(word::AbstractVector{<:Integer}) = join(word, ',')
_word_lt(a::Vector{UInt16}, b::Vector{UInt16}) = isless((length(a), Tuple(a)), (length(b), Tuple(b)))

function _checked_word(word)
    result = UInt16[]
    for label in word
        label isa Integer || throw(ArgumentError("Pauli labels must be integers"))
        1 <= label <= typemax(UInt16) || throw(ArgumentError("Pauli label $label is outside UInt16 site encoding"))
        push!(result, UInt16(label))
    end
    reduced, phase = pauli_product(result)
    return reduced, phase
end

"""Multiply a sequence of encoded Pauli factors without applying physical symmetries."""
function pauli_product(word::AbstractVector{<:Integer}; realify::Bool=false)
    checked = UInt16[]
    for label in word
        label isa Integer || throw(ArgumentError("Pauli labels must be integers"))
        1 <= label <= typemax(UInt16) || throw(ArgumentError("Pauli label $label is outside UInt16 site encoding"))
        push!(checked, UInt16(label))
    end
    isempty(checked) && return UInt16[], 1
    reduce1!(checked)
    reduce3!(checked)
    checked, coefficient = reduce2!(checked, realify=realify)
    reduce3!(checked)
    return checked, coefficient
end

"""A deterministic canonical sum of Pauli terms."""
struct PauliPolynomial
    terms::Vector{PauliTerm}
    hermitian::Bool
end

function PauliPolynomial(terms; hermitian::Bool=true)
    coefficients = Dict{Vector{UInt16},ComplexF64}()
    for item in terms
        word, coefficient = if item isa PauliTerm
            item.word, item.coefficient
        elseif item isa Pair
            item.first, item.second
        elseif item isa Tuple && length(item) == 2
            item[1], item[2]
        else
            throw(ArgumentError("polynomial terms must be PauliTerm, Pair, or (word, coefficient)"))
        end
        isfinite(real(coefficient)) && isfinite(imag(coefficient)) ||
            throw(ArgumentError("Pauli polynomial coefficients must be finite"))
        canonical, phase = _checked_word(word)
        coefficients[canonical] = get(coefficients, canonical, 0.0 + 0.0im) + ComplexF64(coefficient * phase)
    end
    canonical_terms = PauliTerm[]
    for (word, coefficient) in coefficients
        coefficient == 0 && continue
        hermitian && !iszero(imag(coefficient)) &&
            throw(ArgumentError("Hermitian Pauli polynomial has non-real coefficient $coefficient on word $word"))
        push!(canonical_terms, PauliTerm(copy(word), coefficient))
    end
    sort!(canonical_terms, lt=(a, b) -> _word_lt(a.word, b.word))
    return PauliPolynomial(canonical_terms, hermitian)
end
PauliPolynomial() = PauliPolynomial(Pair{Vector{UInt16},Float64}[])

function _poly_dict(polynomial::PauliPolynomial)
    Dict(term.word => term.coefficient for term in polynomial.terms)
end

function _polynomial_product(left::PauliPolynomial, right::PauliPolynomial; hermitian=false)
    terms = Pair{Vector{UInt16},ComplexF64}[]
    for a in left.terms, b in right.terms
        word, phase = pauli_product([a.word; b.word])
        push!(terms, word => a.coefficient * b.coefficient * phase)
    end
    return PauliPolynomial(terms; hermitian=hermitian)
end

function _polynomial_linear_combination(parts; hermitian=false)
    terms = Pair{Vector{UInt16},ComplexF64}[]
    for (scale, polynomial) in parts, term in polynomial.terms
        push!(terms, term.word => scale * term.coefficient)
    end
    return PauliPolynomial(terms; hermitian=hermitian)
end

"""Optional user-declared translation orbit metadata; no orbit is inferred."""
struct TranslationOrbitMetadata
    period::Int
    orbit_ids::Vector{Int}
    function TranslationOrbitMetadata(period, orbit_ids)
        period > 0 || throw(ArgumentError("translation period must be positive"))
        all(>(0), orbit_ids) || throw(ArgumentError("translation orbit ids must be positive"))
        new(period, Int.(orbit_ids))
    end
end

"""A named, explicitly ordered PSD basis sector."""
struct BasisSector
    name::Symbol
    words::Vector{Vector{UInt16}}
    psd_role::Symbol
    translation_orbits::Union{Nothing,TranslationOrbitMetadata}
end

function BasisSector(name, words; psd_role=:moment, translation_orbits=nothing)
    name = Symbol(name)
    isempty(String(name)) && throw(ArgumentError("basis sector name cannot be empty"))
    psd_role in (:moment, :state_optimality, :auxiliary) || throw(ArgumentError("invalid PSD role $psd_role"))
    canonical = Vector{UInt16}[]
    seen = Set{Vector{UInt16}}()
    for word in words
        reduced, phase = _checked_word(word)
        phase == 1 || throw(ArgumentError("basis word $word is not a canonical Hermitian Pauli word"))
        reduced in seen && throw(ArgumentError("duplicate word $reduced in basis sector $name"))
        push!(seen, reduced)
        push!(canonical, reduced)
    end
    isempty(canonical) && throw(ArgumentError("basis sector $name is empty"))
    translation_orbits !== nothing && length(translation_orbits.orbit_ids) != length(canonical) &&
        throw(ArgumentError("translation orbit metadata length does not match sector $name"))
    return BasisSector(name, canonical, psd_role, translation_orbits)
end

"""
An explicit symmetry generator and its single declared use. `site_map[i]` and
`axis_map[a]` define its action; `axis_sign[a]` supplies a ±1 factor.
"""
struct SymmetryDeclaration
    name::Symbol
    purpose::Symbol
    site_map::Vector{Int}
    axis_map::NTuple{3,UInt8}
    axis_sign::NTuple{3,Int8}
end

function SymmetryDeclaration(name, purpose, site_map;
                             axis_map=(1, 2, 3), axis_sign=(1, 1, 1))
    purpose = Symbol(purpose)
    purpose in SYMMETRY_PURPOSES || throw(ArgumentError("invalid symmetry purpose $purpose"))
    sites = Int.(site_map)
    sort(sites) == collect(1:length(sites)) || throw(ArgumentError("site_map must be a permutation of 1:L"))
    amap = Tuple(UInt8.(axis_map))
    sort(collect(amap)) == UInt8[1, 2, 3] || throw(ArgumentError("axis_map must permute X,Y,Z"))
    signs = Tuple(Int8.(axis_sign))
    all(sign -> sign in (-1, 1), signs) || throw(ArgumentError("axis signs must be ±1"))
    return SymmetryDeclaration(Symbol(name), purpose, sites, amap, signs)
end

function _transform_word(word::Vector{UInt16}, symmetry::SymmetryDeclaration)
    transformed = UInt16[]
    coefficient = 1
    for label in word
        site = cld(Int(label), 3)
        site <= length(symmetry.site_map) || throw(ArgumentError("symmetry $(symmetry.name) does not cover site $site"))
        axis = Int(mod1(label, 3))
        coefficient *= symmetry.axis_sign[axis]
        push!(transformed, UInt16(3 * (symmetry.site_map[site] - 1) + symmetry.axis_map[axis]))
    end
    canonical, phase = pauli_product(transformed)
    return canonical, coefficient * phase
end

function _check_invariance(polynomial::PauliPolynomial, symmetry::SymmetryDeclaration)
    transformed = Pair{Vector{UInt16},ComplexF64}[]
    for term in polynomial.terms
        word, coefficient = _transform_word(term.word, symmetry)
        push!(transformed, word => coefficient * term.coefficient)
    end
    transformed_polynomial = PauliPolynomial(transformed; hermitian=polynomial.hermitian)
    _poly_dict(transformed_polynomial) == _poly_dict(polynomial) ||
        throw(ArgumentError("Hamiltonian is not invariant under declared symmetry $(symmetry.name) ($(symmetry.purpose))"))
end

"""An explicitly ordered subsystem for a mechanically expanded small RDM."""
struct RDMRegion
    name::Symbol
    sites::Vector{Int}
    blocks::Union{Nothing,Vector{Vector{Int}}}
end
function RDMRegion(name, sites; blocks=nothing)
    sites = Int.(sites)
    !isempty(sites) && all(>(0), sites) && length(unique(sites)) == length(sites) ||
        throw(ArgumentError("RDM sites must be distinct positive integers"))
    length(sites) <= 4 || throw(ArgumentError("Gate A explicit RDM expansion supports at most four sites"))
    return RDMRegion(Symbol(name), sites, blocks)
end

function _maximum_site(polynomial::PauliPolynomial)
    maximum((cld(Int(label), 3) for term in polynomial.terms for label in term.word); init=0)
end

"""Complete explicit input to the Gate A compiler."""
struct RelaxationSpecification
    hamiltonian::PauliPolynomial
    basis::Vector{BasisSector}
    symmetries::Vector{SymmetryDeclaration}
    linear_tests::Vector{PauliPolynomial}
    psd_state_basis::Vector{Vector{UInt16}}
    rdm_regions::Vector{RDMRegion}
    observables::Dict{Symbol,PauliPolynomial}
    normalization::Float64
    certificate_scope::Symbol
    declared_moment_support::Union{Nothing,Vector{Vector{UInt16}}}
end

function RelaxationSpecification(hamiltonian, basis;
        symmetries=SymmetryDeclaration[], linear_tests=PauliPolynomial[],
        psd_state_basis=Vector{UInt16}[], rdm_regions=RDMRegion[],
        observables=Dict{Symbol,PauliPolynomial}(), normalization=1.0,
        certificate_scope=:numerical_relaxation, declared_moment_support=nothing)
    hamiltonian isa PauliPolynomial && hamiltonian.hermitian || throw(ArgumentError("Hamiltonian must be a Hermitian PauliPolynomial"))
    isempty(basis) && throw(ArgumentError("at least one explicit basis sector is required"))
    names = getfield.(basis, :name)
    length(unique(names)) == length(names) || throw(ArgumentError("basis sector names must be unique"))
    all(test -> test isa PauliPolynomial, linear_tests) ||
        throw(ArgumentError("linear tests must be explicit PauliPolynomial values"))
    all(observable -> observable isa PauliPolynomial && observable.hermitian, values(observables)) ||
        throw(ArgumentError("observables must be Hermitian PauliPolynomial values"))
    all(region -> region isa RDMRegion, rdm_regions) ||
        throw(ArgumentError("RDM regions must be explicit RDMRegion values"))
    scope = validate_relaxation_label(:result_scope, certificate_scope)
    scope == :rigorously_postvalidated && throw(ArgumentError("rigorously_postvalidated scope requires an independent postvalidation result"))
    isfinite(normalization) && normalization != 0 || throw(ArgumentError("normalization must be finite and nonzero"))
    maxsite = _maximum_site(hamiltonian)
    for sector in basis, word in sector.words
        for label in word
            maxsite = max(maxsite, cld(Int(label), 3))
        end
    end
    for test in linear_tests
        maxsite = max(maxsite, _maximum_site(test))
    end
    for observable in values(observables)
        maxsite = max(maxsite, _maximum_site(observable))
    end
    for word in psd_state_basis, label in word
        maxsite = max(maxsite, cld(Int(label), 3))
    end
    for region in rdm_regions
        maxsite = max(maxsite, maximum(region.sites))
    end
    for symmetry in symmetries
        length(symmetry.site_map) >= maxsite || throw(ArgumentError("symmetry $(symmetry.name) does not cover all used sites"))
        _check_invariance(hamiltonian, symmetry)
    end
    psd_words = Vector{UInt16}[]
    for word in psd_state_basis
        canonical, phase = _checked_word(word)
        phase == 1 || throw(ArgumentError("PSD state-optimality basis contains a phased word"))
        push!(psd_words, canonical)
    end
    support = declared_moment_support === nothing ? nothing : [_checked_word(word)[1] for word in declared_moment_support]
    return RelaxationSpecification(hamiltonian, collect(basis), collect(symmetries), collect(linear_tests),
        psd_words, collect(rdm_regions), Dict{Symbol,PauliPolynomial}(observables), Float64(normalization), scope, support)
end

struct AffineMomentEntry
    words::Vector{Vector{UInt16}}
    coefficients::Vector{ComplexF64}
end

struct CompiledPSDBlock
    name::Symbol
    role::Symbol
    entries::Matrix{AffineMomentEntry}
    provenance::Vector{String}
end

struct CompilationDiagnostics
    raw_basis_size::Int
    canonical_basis_size::Int
    raw_bdagb_entries::Int
    zeroed_moments::Int
    equated_moments::Int
    scalar_moment_count::Int
    linear_rows_raw::Int
    linear_rows_zero::Int
    linear_rows_duplicate::Int
    linear_rows_independent::Int
    psd_block_count::Int
    psd_block_dimensions::Vector{Int}
    max_psd_block_dimension::Int
    real_embedded_dimensions::Vector{Int}
    rdm_moment_increment::Int
    state_optimality_moment_increment::Int
    missing_support::Vector{String}
    constraint_provenance::Vector{String}
    scope::Symbol
    fingerprint::String
end

struct CompiledRelaxation
    specification::RelaxationSpecification
    moments::Vector{Vector{UInt16}}
    moment_index::Dict{Vector{UInt16},Int}
    objective::AffineMomentEntry
    observables::Dict{Symbol,AffineMomentEntry}
    psd_blocks::Vector{CompiledPSDBlock}
    linear_rows::Vector{AffineMomentEntry}
    diagnostics::CompilationDiagnostics
end

function _canonical_moment(word::Vector{UInt16}, symmetries)
    generators = filter(symmetry -> symmetry.purpose in (:moment_zero, :moment_equality), symmetries)
    isempty(generators) && return copy(word), 1.0 + 0.0im, :unchanged

    phases = Dict{Vector{UInt16},ComplexF64}(copy(word) => 1.0 + 0.0im)
    queue = Vector{Vector{UInt16}}([copy(word)])
    cursor = 1
    while cursor <= length(queue)
        current = queue[cursor]
        current_phase = phases[current]
        cursor += 1
        for symmetry in generators
            transformed, factor = _transform_word(current, symmetry)
            transformed_phase = current_phase * factor
            if haskey(phases, transformed)
                phases[transformed] == transformed_phase ||
                    return UInt16[], 0.0 + 0.0im, :symmetry_zero
            else
                phases[transformed] = transformed_phase
                push!(queue, transformed)
            end
        end
    end
    canonical = first(sort!(collect(keys(phases)); lt=_word_lt))
    return copy(canonical), phases[canonical], canonical == word ? :unchanged : :equated
end

function _entry(polynomial::PauliPolynomial, symmetries; provenance="")
    combined = Dict{Vector{UInt16},ComplexF64}()
    zeroed = equated = 0
    for term in polynomial.terms
        word, factor, reason = _canonical_moment(term.word, symmetries)
        if factor == 0
            zeroed += 1
        else
            reason == :equated && (equated += 1)
            combined[word] = get(combined, word, 0.0 + 0.0im) + factor * term.coefficient
        end
    end
    words = sort!(collect(keys(combined)), lt=_word_lt)
    filter!(word -> combined[word] != 0, words)
    return AffineMomentEntry(words, ComplexF64[combined[word] for word in words]), zeroed, equated
end

_entry_for_product(left, right, symmetries) = _entry(_polynomial_product(
    PauliPolynomial([left => 1.0]), PauliPolynomial([right => 1.0]); hermitian=false), symmetries)

function _row_signature(entry::AffineMomentEntry)
    isempty(entry.words) && return ""
    first_nonzero = findfirst(!iszero, entry.coefficients)
    first_nonzero === nothing && return ""
    scale = entry.coefficients[first_nonzero]
    join((_word_key(word) * "=" * repr(coefficient / scale)
          for (word, coefficient) in zip(entry.words, entry.coefficients)), ';')
end

function _pauli_matrix(axis::Int)
    axis == 0 && return ComplexF64[1 0; 0 1]
    axis == 1 && return ComplexF64[0 1; 1 0]
    axis == 2 && return ComplexF64[0 -im; im 0]
    return ComplexF64[1 0; 0 -1]
end

function _rdm_entries(region::RDMRegion, symmetries)
    k = length(region.sites)
    dimension = 2^k
    accum = [Dict{Vector{UInt16},ComplexF64}() for _ in 1:dimension, _ in 1:dimension]
    zeroed = equated = 0
    for axes in Iterators.product(ntuple(_ -> 0:3, k)...)
        word = UInt16[3 * (region.sites[i] - 1) + axes[i] for i in 1:k if axes[i] != 0]
        canonical, factor, reason = _canonical_moment(word, symmetries)
        factor == 0 && (zeroed += 1; continue)
        reason == :equated && (equated += 1)
        matrix = _pauli_matrix(axes[1])
        for i in 2:k
            matrix = kron(matrix, _pauli_matrix(axes[i]))
        end
        matrix ./= 2^k
        for row in 1:dimension, column in 1:dimension
            coefficient = factor * matrix[row, column]
            coefficient == 0 && continue
            accum[row, column][canonical] = get(accum[row, column], canonical, 0.0 + 0.0im) + coefficient
        end
    end
    entries = Matrix{AffineMomentEntry}(undef, dimension, dimension)
    for row in 1:dimension, column in 1:dimension
        words = sort!(collect(keys(accum[row, column])), lt=_word_lt)
        entries[row, column] = AffineMomentEntry(words, [accum[row, column][word] for word in words])
    end
    return entries, zeroed, equated
end

function _fingerprint(text::AbstractString)
    hash = UInt64(0xcbf29ce484222325)
    for byte in codeunits(text)
        hash = (hash ⊻ UInt64(byte)) * UInt64(0x100000001b3)
    end
    return lowercase(string(hash, base=16, pad=16))
end

function _artifact_fingerprint(spec, moments, blocks, rows, objective)
    io = IOBuffer()
    print(io, "scope=", spec.certificate_scope, ";normalization=", repr(spec.normalization), ';')
    for term in spec.hamiltonian.terms
        print(io, "H:", _word_key(term.word), '=', repr(term.coefficient), ';')
    end
    for word in moments
        print(io, "M:", _word_key(word), ';')
    end
    for block in blocks
        print(io, "P:", block.name, ':', size(block.entries), ';')
        for entry in block.entries, (word, coefficient) in zip(entry.words, entry.coefficients)
            print(io, _word_key(word), '=', repr(coefficient), ';')
        end
    end
    for row in rows
        print(io, "L:", _row_signature(row), ';')
    end
    print(io, "O:", _row_signature(objective))
    return _fingerprint(String(take!(io)))
end

function _scaled_entry_sum(parts)
    combined = Dict{Vector{UInt16},ComplexF64}()
    for (scale, entry) in parts, (word, coefficient) in zip(entry.words, entry.coefficients)
        combined[word] = get(combined, word, 0.0 + 0.0im) + scale * coefficient
    end
    words = sort!(collect(keys(combined)), lt=_word_lt)
    filter!(word -> combined[word] != 0, words)
    return AffineMomentEntry(words, ComplexF64[combined[word] for word in words])
end

function _basis_block_indices(sector::BasisSector, declarations)
    generators = filter(declaration -> declaration.purpose == :basis_block, declarations)
    isempty(generators) && return [(sector.name, collect(eachindex(sector.words)))]
    groups = Dict{Tuple,Vector{Int}}()
    for (index, word) in enumerate(sector.words)
        signature = Tuple(begin
            transformed, factor = _transform_word(word, generator)
            transformed == word || throw(ArgumentError(
                "basis_block generator $(generator.name) does not act diagonally on $(sector.name) word $index"))
            isreal(factor) && real(factor) in (-1, 1) || throw(ArgumentError(
                "basis_block generator $(generator.name) has non-character phase on $(sector.name) word $index"))
            Int(real(factor))
        end for generator in generators)
        push!(get!(groups, signature, Int[]), index)
    end
    signatures = sort!(collect(keys(groups)))
    return [(Symbol(sector.name, :_character_, join(signature, '_')), groups[signature])
            for signature in signatures]
end

function _fourier_blocks(sector::BasisSector, entries, declarations)
    generators = filter(declaration -> declaration.purpose == :fourier_orbit, declarations)
    isempty(generators) && return nothing
    length(generators) == 1 || throw(ArgumentError("Gate A supports one explicit Fourier generator per sector"))
    metadata = sector.translation_orbits
    metadata === nothing && throw(ArgumentError("Fourier sector $(sector.name) requires explicit translation orbit metadata"))
    period = metadata.period
    generator = only(generators)
    orbit_ids = unique(metadata.orbit_ids)
    orbits = [findall(==(orbit_id), metadata.orbit_ids) for orbit_id in orbit_ids]
    all(length(orbit) == period for orbit in orbits) ||
        throw(ArgumentError("Fourier sector $(sector.name) must list complete orbits of period $period"))
    for orbit in orbits, position in 1:period
        transformed, factor = _transform_word(sector.words[orbit[position]], generator)
        expected = sector.words[orbit[mod1(position + 1, period)]]
        transformed == expected && factor == 1 || throw(ArgumentError(
            "Fourier generator $(generator.name) does not advance declared orbit in $(sector.name)"))
    end
    blocks = CompiledPSDBlock[]
    for momentum in 0:(period - 1)
        dimension = length(orbits)
        transformed_entries = Matrix{AffineMomentEntry}(undef, dimension, dimension)
        for left in 1:dimension, right in 1:dimension
            parts = Tuple{ComplexF64,AffineMomentEntry}[]
            for left_position in 1:period, right_position in 1:period
                phase = cis(2pi * momentum * (left_position - right_position) / period) / period
                push!(parts, (phase, entries[orbits[left][left_position], orbits[right][right_position]]))
            end
            transformed_entries[left, right] = _scaled_entry_sum(parts)
        end
        push!(blocks, CompiledPSDBlock(Symbol(sector.name, :_momentum_, momentum),
            sector.psd_role, transformed_entries,
            ["Fourier $(generator.name), momentum $momentum [$row,$column]"
             for row in 1:dimension for column in 1:dimension]))
    end
    return blocks
end

"""Compile an explicit relaxation without constructing or invoking an optimizer."""
function compile_relaxation(spec::RelaxationSpecification)
    blocks = CompiledPSDBlock[]
    all_words = Set{Vector{UInt16}}([UInt16[]])
    zeroed = equated = 0
    provenance = String[]
    raw_basis_size = sum(length(sector.words) for sector in spec.basis)

    for sector in spec.basis
        n = length(sector.words)
        entries = Matrix{AffineMomentEntry}(undef, n, n)
        for row in 1:n, column in 1:n
            product_word, phase = pauli_product([reverse(sector.words[row]); sector.words[column]])
            polynomial = PauliPolynomial([product_word => phase]; hermitian=false)
            entry, z, e = _entry(polynomial, spec.symmetries; provenance="B†B $(sector.name)[$row,$column]")
            entries[row, column] = entry
            union!(all_words, entry.words)
            zeroed += z; equated += e
        end
        fourier_blocks = _fourier_blocks(sector, entries, spec.symmetries)
        if fourier_blocks !== nothing
            append!(blocks, fourier_blocks)
            append!(provenance, ["Fourier PSD sector $(block.name), dimension $(size(block.entries, 1))"
                                 for block in fourier_blocks])
        else
            for (block_name, indices) in _basis_block_indices(sector, spec.symmetries)
                block_entries = entries[indices, indices]
                push!(blocks, CompiledPSDBlock(block_name, sector.psd_role, block_entries,
                    ["B†B:$(sector.name)[$row,$column]" for row in indices for column in indices]))
                push!(provenance, "moment PSD sector $block_name, dimension $(length(indices))")
            end
        end
    end

    objective, z, e = _entry(spec.hamiltonian, spec.symmetries; provenance="Hamiltonian objective")
    zeroed += z; equated += e; union!(all_words, objective.words)
    observables = Dict{Symbol,AffineMomentEntry}()
    for name in sort!(collect(keys(spec.observables)))
        entry, z, e = _entry(spec.observables[name], spec.symmetries; provenance="observable $name")
        observables[name] = entry; zeroed += z; equated += e; union!(all_words, entry.words)
    end

    raw_linear = length(spec.linear_tests)
    zero_linear = duplicate_linear = 0
    linear_rows = AffineMomentEntry[]
    signatures = Set{String}()
    for (test_index, test) in enumerate(spec.linear_tests)
        commutator = _polynomial_linear_combination(((1.0, _polynomial_product(spec.hamiltonian, test)),
                                                     (-1.0, _polynomial_product(test, spec.hamiltonian))); hermitian=false)
        entry, z, e = _entry(commutator, spec.symmetries; provenance="linear test $test_index")
        zeroed += z; equated += e
        signature = _row_signature(entry)
        if isempty(signature)
            zero_linear += 1
        elseif signature in signatures
            duplicate_linear += 1
        else
            push!(signatures, signature); push!(linear_rows, entry); union!(all_words, entry.words)
            push!(provenance, "linear state optimality test $test_index")
        end
    end

    before_pso = length(all_words)
    if !isempty(spec.psd_state_basis)
        n = length(spec.psd_state_basis)
        entries = Matrix{AffineMomentEntry}(undef, n, n)
        for row in 1:n, column in row:n
            v = PauliPolynomial([spec.psd_state_basis[row] => 1.0])
            wdag = PauliPolynomial([reverse(spec.psd_state_basis[column]) => 1.0])
            vwdag = _polynomial_product(v, wdag)
            expression = _polynomial_linear_combination((
                (1.0, _polynomial_product(_polynomial_product(v, spec.hamiltonian), wdag)),
                (-0.5, _polynomial_product(spec.hamiltonian, vwdag)),
                (-0.5, _polynomial_product(vwdag, spec.hamiltonian))), hermitian=false)
            entry, z, e = _entry(expression, spec.symmetries; provenance="PSD state optimality[$row,$column]")
            entries[row, column] = entry
            if row != column
                entries[column, row] = AffineMomentEntry(copy(entry.words), conj.(entry.coefficients))
            end
            zeroed += z; equated += e; union!(all_words, entry.words)
        end
        push!(blocks, CompiledPSDBlock(:state_optimality, :state_optimality, entries,
            ["PSD state optimality[$row,$column]" for row in 1:n for column in 1:n]))
        push!(provenance, "PSD state optimality, dimension $n")
    end
    pso_increment = length(all_words) - before_pso

    before_rdm = length(all_words)
    for region in spec.rdm_regions
        entries, z, e = _rdm_entries(region, spec.symmetries)
        zeroed += z; equated += e
        for entry in entries; union!(all_words, entry.words); end
        declared_blocks = region.blocks === nothing ? [collect(axes(entries, 1))] : region.blocks
        flattened = vcat(declared_blocks...)
        sort(flattened) == collect(axes(entries, 1)) && length(unique(flattened)) == length(flattened) ||
            throw(ArgumentError("RDM $(region.name) blocks must partition matrix indices 1:$(size(entries, 1))"))
        for (block_index, indices) in enumerate(declared_blocks)
            block_entries = entries[indices, indices]
            block_name = length(declared_blocks) == 1 ? region.name : Symbol(region.name, :_block_, block_index)
            push!(blocks, CompiledPSDBlock(block_name, :rdm, block_entries,
                ["RDM $(region.name) block $block_index [$row,$column]" for row in indices for column in indices]))
            push!(provenance, "RDM $(region.name) block $block_index, dimension $(length(indices))")
        end
    end
    rdm_increment = length(all_words) - before_rdm

    moments = sort!(collect(all_words), lt=_word_lt)
    missing = String[]
    if spec.declared_moment_support !== nothing
        declared = Set(spec.declared_moment_support)
        for word in moments
            word in declared || push!(missing, "missing canonical moment [$(_word_key(word))]")
        end
        isempty(missing) || throw(ArgumentError("moment closure failed: " * join(missing, "; ")))
    end
    index = Dict(word => i for (i, word) in enumerate(moments))

    # Deterministic algebraic Hermiticity checks for every PSD artifact.
    for block in blocks, row in axes(block.entries, 1), column in axes(block.entries, 2)
        lhs = Dict(zip(block.entries[row, column].words, block.entries[row, column].coefficients))
        rhs = Dict(word => conj(coefficient) for (word, coefficient) in
                   zip(block.entries[column, row].words, block.entries[column, row].coefficients))
        lhs == rhs || throw(ArgumentError("PSD Hermiticity failed in $(block.name) at ($row,$column)"))
    end
    isempty(blocks) && throw(ArgumentError("relaxation has no PSD block"))
    any(block -> size(block.entries, 1) == 0, blocks) && throw(ArgumentError("empty PSD block"))

    dimensions = [size(block.entries, 1) for block in blocks]
    fingerprint = _artifact_fingerprint(spec, moments, blocks, linear_rows, objective)
    diagnostics = CompilationDiagnostics(raw_basis_size, raw_basis_size,
        sum(length(sector.words)^2 for sector in spec.basis), zeroed, equated,
        length(moments), raw_linear, zero_linear, duplicate_linear, length(linear_rows),
        length(blocks), dimensions, maximum(dimensions), 2 .* dimensions,
        rdm_increment, pso_increment, missing, provenance, spec.certificate_scope, fingerprint)
    return CompiledRelaxation(spec, moments, index, objective, observables, blocks, linear_rows, diagnostics)
end

struct BuiltRelaxationModel
    model::Model
    moments
    compiled::CompiledRelaxation
end

function _jump_entry(entry::AffineMomentEntry, moments)
    real_expression = AffExpr(0.0)
    imaginary_expression = AffExpr(0.0)
    for (word, coefficient) in zip(entry.words, entry.coefficients)
        index = moments.index[word]
        add_to_expression!(real_expression, real(coefficient), moments.variables[index])
        add_to_expression!(imaginary_expression, imag(coefficient), moments.variables[index])
    end
    return real_expression, imaginary_expression
end

"""Build JuMP variables and constraints from a compiled artifact, without solving it."""
function build_jump_model(compiled::CompiledRelaxation; optimizer=nothing, optimizer_attributes=Pair[])
    model = optimizer === nothing ? Model() : Model(optimizer_with_attributes(optimizer, optimizer_attributes...))
    @variable(model, moment_variables[1:length(compiled.moments)])
    moment_handle = (variables=moment_variables, index=compiled.moment_index)
    @constraint(model, moment_variables[compiled.moment_index[UInt16[]]] == 1)
    for row in compiled.linear_rows
        real_row, imaginary_row = _jump_entry(row, moment_handle)
        @constraint(model, real_row == 0)
        @constraint(model, imaginary_row == 0)
    end
    for block in compiled.psd_blocks
        n = size(block.entries, 1)
        real_part = Matrix{AffExpr}(undef, n, n)
        imaginary_part = Matrix{AffExpr}(undef, n, n)
        for row in 1:n, column in 1:n
            real_part[row, column], imaginary_part[row, column] = _jump_entry(block.entries[row, column], moment_handle)
        end
        @constraint(model, Symmetric([real_part -imaginary_part; imaginary_part real_part]) in PSDCone())
    end
    objective, imaginary_objective = _jump_entry(compiled.objective, moment_handle)
    isempty(imaginary_objective.terms) && imaginary_objective.constant == 0 ||
        throw(ArgumentError("Hamiltonian objective is not real"))
    @objective(model, Min, objective / compiled.specification.normalization)
    return BuiltRelaxationModel(model, moment_handle, compiled)
end

struct RelaxationSolveResult
    objective_value::Float64
    status
    primal_status
    scope::Symbol
    fingerprint::String
end

"""Optimize a previously built model and report an explicitly scoped numerical result."""
function solve_relaxation(built::BuiltRelaxationModel)
    optimize!(built.model)
    status = termination_status(built.model)
    has_values(built.model) || error("relaxation terminated without a solution: $status")
    scope = built.compiled.specification.certificate_scope == :numerical_relaxation ? :solver_bound : built.compiled.specification.certificate_scope
    return RelaxationSolveResult(objective_value(built.model), status, primal_status(built.model),
                                 scope, built.compiled.diagnostics.fingerprint)
end

"""Construct the paper's explicit 1D sparse basis and Table 2 structural counts."""
function heisenberg_table2_benchmark(N::Int=100, d::Int=4, r::Int=1)
    N > 0 && 0 <= d <= N && r >= 1 || throw(ArgumentError("invalid Table 2 parameters"))
    words = Vector{UInt16}[UInt16[]]
    for length_ in 1:d, start in 1:N
        for axes in Iterators.product(ntuple(_ -> 1:3, length_)...)
            push!(words, UInt16[3 * (mod(start + offset - 2, N)) + axes[offset] for offset in 1:length_])
        end
    end
    if d >= 2 && r >= 2
        for distance in 2:r, start in 1:N, left_axis in 1:3, right_axis in 1:3
            push!(words, UInt16[3 * (start - 1) + left_axis,
                                3 * mod(start + distance - 1, N) + right_axis])
        end
    end
    basis = BasisSector(:table2_sparse_basis, words; psd_role=:moment)
    sparse_size = length(words)
    max_block = iseven(d) ? div(3^(d + 1) + 5, 8) : div(3^(d + 1) - 1, 8)
    return (basis=basis, sparse_basis_size=sparse_size,
            max_psd_block_dimension=max_block,
            original_dimension=div((3N)^(d + 1) - 1, 3N - 1),
            equality_reduced_dimension=sum(binomial(N, i) * 3^i for i in 0:d))
end

function add_SU2_equality!(model, tsupp, cons; L=0, lattice="chain")
    ind = findall(item->length(item) == 4 && all(smod.(item, 3) .== 1), tsupp)
    for item in tsupp[ind]
        fr = @variable(model)
        Locb = bfind(tsupp, item)
        add_to_expression!(cons[Locb], fr)
        for i = 2:4
            a = copy(item)
            a[1] += 1
            a[i] += 1
            a = reduce!(a, L=L, lattice=lattice)[1]
            Locb = bfind(tsupp, a)
            add_to_expression!(cons[Locb], -1, fr)
        end
    end
    ind = findall(item->length(item) == 6 && sum(smod.(item, 3) .== 1) == 4, tsupp)
    for item in tsupp[ind]
        ino = Vector(1:6)[smod.(item, 3) .== 1]
        fr = @variable(model)
        Locb = bfind(tsupp, item)
        add_to_expression!(cons[Locb], fr)
        for i = 2:4
            a = copy(item)
            a[ino[1]] += 2
            a[ino[i]] += 2
            a = reduce!(a, L=L, lattice=lattice)[1]
            Locb = bfind(tsupp, a)
            add_to_expression!(cons[Locb], -1, fr)
        end
    end
    ind = findall(item->length(item) == 6 && all(smod.(item, 3) .== 1), tsupp)
    for item in tsupp[ind]
        fr = @variable(model)
        Locb = bfind(tsupp, item)
        add_to_expression!(cons[Locb], fr)
        for i = 2:6
            ino = [Vector(2:i-1); Vector(i+1:6)]
            for j = 2:4
                ine = [Vector(2:j-1); Vector(j+1:4)]
                a = copy(item)
                a[ino[1]] += 1
                a[ino[j]] += 1
                a[ino[ine[1]]] += 2
                a[ino[ine[2]]] += 2
                a = reduce!(a, L=L, lattice=lattice)[1]
                Locb = bfind(tsupp, a)
                add_to_expression!(cons[Locb], -1, fr)
            end
        end
    end
    ind = findall(item->length(item) == 8 && sum(smod.(item, 3) .== 1) == 6, tsupp)
    for item in tsupp[ind]
        ino = Vector(1:8)[smod.(item, 3) .== 1]
        for i = 1:6
            fr = @variable(model)
            Locb = bfind(tsupp, item)
            add_to_expression!(cons[Locb], fr)
            for j in [Vector(1:i-1); Vector(i+1:6)]
                a = copy(item)
                a[ino[i]] += 2
                a[ino[j]] += 2
                a = reduce!(a, L=L, lattice=lattice)[1]
                Locb = bfind(tsupp, a)
                add_to_expression!(cons[Locb], -1, fr)
            end
        end
    end
    ind = findall(item->length(item) == 8 && sum(smod.(item, 3) .== 1) == 4 && sum(smod.(item, 3) .== 2) == 4, tsupp)
    for item in tsupp[ind]
        ino = Vector(1:8)[smod.(item, 3) .== 2]
        fr = @variable(model)
        Locb = bfind(tsupp, item)
        add_to_expression!(cons[Locb], fr)
        for i in 2:4
            a = copy(item)
            a[ino[1]] += 1
            a[ino[i]] += 1
            a = reduce!(a, L=L, lattice=lattice)[1]
            Locb = bfind(tsupp, a)
            add_to_expression!(cons[Locb], -1, fr)
        end
    end
end
