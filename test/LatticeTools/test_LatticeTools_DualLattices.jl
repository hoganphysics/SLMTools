using Test
using SLMTools.LatticeTools.DualLattices
# using SLMTools

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


# Define sample lattices and LatticeFields for testing
L1 = (1:10, 1:20)
L2 = dualShiftLattice(L1)  # Creating a dual lattice for comparison
flambda = 1.0

data = rand(ComplexF64, 10, 20)
data2 = rand(ComplexF64, 10, 20)
LF1 = LF{ComplexPhase}(data, L1, flambda)
LF2 = LF{ComplexPhase}(data, L2, flambda)
LF3 = LF{ComplexPhase}(data2, L2, flambda)

@testset "Lattice Duality Query Tests" begin
    # Test ldq for lattices
    @testset "ldq for Lattices" begin
        @test ldq(L1, L2, flambda) == nothing
        @test_throws DimensionMismatch ldq(L1, (1:10, 1:10), flambda)
        @test_throws DomainError ldq(L1, L1, flambda)
    end

    # Test ldq for LatticeFields
    @testset "ldq for LatticeFields" begin
        @test ldq(LF1, LF2) == nothing
        @test_throws DimensionMismatch ldq(LF1, LF{ComplexPhase}(data, (1:10, 1:10), flambda))
        @test_throws DomainError ldq(LF1, LF1)
    end
end



@testset "Dual Phase Tests" begin
    @testset "dualPhase without specifying dL" begin
        dpField = dualPhase(L1)
        @test isa(dpField, LF{ComplexPhase})
        @test dpField.L == dualShiftLattice(L1)
    end

    @testset "dualPhase with specifying dL" begin
        dL = dualShiftLattice(L1)
        dpField = dualPhase(L1, dL)
        @test isa(dpField, LF{ComplexPhase})
        @test dpField.L == dL
    end
end
