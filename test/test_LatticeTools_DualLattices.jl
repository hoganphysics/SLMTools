using Test
using SLMTools.LatticeTools.DualLattices

@testset "DualLattices tests" begin
    @testset "dualLattice" begin
        L = (1:10, 1:20)
        flambda = 2.0
        dL = dualLattice(L, flambda)
        isapprox_range1 = all(isapprox.(dL[1], (0.0:0.2:1.8), atol=1e-10))
        isapprox_range2 = all(isapprox.(dL[2], (0.0:0.1:1.9000000000000001), atol=1e-10))
        @test result = isapprox_range1 && isapprox_range2
    end

    @testset "dualShiftLattice" begin
        L = (1:10, 1:20)
        flambda = 2.0
        dL = dualShiftLattice(L, flambda)
        isapprox_range1 = all(isapprox.(dL[1], (-1.0:0.2:0.8), atol=1e-10))
        isapprox_range2 = all(isapprox.(dL[2], (-1.0:0.1:0.9), atol=1e-10))
        @test result = isapprox_range1 && isapprox_range2
    end
end