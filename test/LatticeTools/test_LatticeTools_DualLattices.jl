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

# Define a sample lattice for testing
sampleLattice = (1:10, 1:20)
flambda = 2.0

@testset "Dual Lattice Tests" begin
    # Test for dualLattice function
    @testset "Testing dualLattice" begin
        dualLat = dualLattice(sampleLattice, flambda)
        @test length(dualLat) == length(sampleLattice)

        # Check step sizes for each dimension
        for (l, dl) in zip(sampleLattice, dualLat)
            expected_step_size = flambda / (length(l) * step(l))
            actual_step_size = step(dl)
            @test actual_step_size ≈ expected_step_size
        end
    end

    # Test for dualShiftLattice function
    @testset "Testing dualShiftLattice" begin
        shiftDualLat = dualShiftLattice(sampleLattice, flambda)
        @test length(shiftDualLat) == length(sampleLattice)

        # Check step sizes and range starts for each dimension
        for (l, sdl) in zip(sampleLattice, shiftDualLat)
            expected_step_size = flambda / (length(l) * step(l))
            actual_step_size = step(sdl)
            @test actual_step_size ≈ expected_step_size

            # Check if the initial entry is negative (for non-even lengths)
            if length(l) % 2 != 0
                @test first(sdl) < 0
            end
        end
    end
end
