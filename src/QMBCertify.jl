module QMBCertify

import Base: iszero

using JuMP
using MosekTools
using LinearAlgebra
using Dualization
using DynamicPolynomials

include(joinpath(@__DIR__, "certification", "helpers.jl"))
include(joinpath(@__DIR__, "certification", "energy_cert.jl"))
include(joinpath(@__DIR__, "certification", "corr_cert.jl"))

export certify_qmb, certify_qmb_corr, dmrg_heisenberg_rat

export GSB, PFB, slabel, reduce!, mosek_para
export AbstractSymmetryModel, PauliSymmetryModel, heisenberg_symmetry, ising_chain_symmetry
export symmetry_reduce_support
export RELAXATION_TERMINOLOGY, validate_relaxation_label
export PauliTerm, PauliPolynomial, pauli_product
export TranslationOrbitMetadata, BasisSector, SymmetryDeclaration, RDMRegion
export RelaxationSpecification, AffineMomentEntry, CompiledPSDBlock
export CompilationDiagnostics, CompiledRelaxation, BuiltRelaxationModel, RelaxationSolveResult
export compile_relaxation, build_jump_model, solve_relaxation, heisenberg_table2_benchmark
export dimerized_j1j2_hamiltonian, dimerized_chain_basis, mg_three_site_projector
export dimerized_chain_observables, dimerized_chain_exact_benchmark, dimerized_chain_scan
export dimerized_chain_symmetries, dimerized_chain_specification
export cluster_chain_hamiltonian, cluster_chain_basis, cluster_chain_observables
export cluster_chain_symmetries, cluster_chain_specification
export IsingMomentResult, ising_relaxation_specification, compile_ising_relaxation
export ising_ground_state_bound

mutable struct qmb_data
    correlation1
    correlation2
    correlation3
    basis
    sbasis
    tsupp
    GramMat
    sGramMat
    multiplier
    moment
end

include("basic_function.jl")
include("ising.jl")
include("rdm_positivity.jl")
include("bound_gsp.jl")
include("bound_partfunc.jl")



end
