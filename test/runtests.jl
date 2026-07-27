using Test
using QMBCertify

@testset "model-specific Pauli symmetry" begin
    heisenberg = heisenberg_symmetry()
    ising = ising_chain_symmetry()

    @test reduce!(UInt16[1], L=4, symmetry=heisenberg)[2] == 0
    @test reduce!(UInt16[1], L=4, symmetry=ising) == (UInt16[1], 1)
    @test reduce!(UInt16[3], L=4, symmetry=ising)[2] == 0
    @test reduce!(UInt16[3, 6], L=4, symmetry=ising)[2] != 0
    @test reduce!(UInt16[4], L=4, symmetry=ising)[1] == UInt16[1]

    reduced = symmetry_reduce_support([[1], [4], [3], [3, 6]], 4, symmetry=ising)
    @test reduced == [UInt16[1], UInt16[3, 6]]
end

@testset "exact transverse-field Ising endpoints" begin
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
end
