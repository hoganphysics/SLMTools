using Test
using SLMTools.LatticeTools

include("test_LatticeTools_Core.jl")
include("test_LatticeTools_DualLattices.jl")
include("test_LatticeTools_Resampling.jl")
# @testset "dualLattice and dualShiftLattice" begin
#     # Test 1: Check that dualLattice works correctly
#     L = (1:5, 1:5)
#     flambda = 2
#     result = dualLattice(L, flambda)
#     @test size(result[1]) == (5,)
#     @test size(result[2]) == (5,)

#     # Test 2: Check that dualShiftLattice works correctly
#     result = dualShiftLattice(L, flambda)
#     @test size(result[1]) == (5,)
#     @test size(result[2]) == (5,)

#     # Test 3: Check that padout works correctly
#     L = 1:5
#     p = (2, 3)
#     result = padout(L, p)
#     @test length(result) == 10
#     @test result[1] == -1
#     @test result[end] == 8
# end
# @testset "dualate function" begin
#     # Test 1: Check that the function throws an error when the camera grid and image size are unequal
#     @test_throws ErrorException dualate([rand(5, 5)], (1:5, 1:5), (1:4, 1:4), 0, 1.0)

#     # Test 2: Check that the function works correctly with a single image and a square grid
#     img = [rand(5, 5)]
#     slmGrid = (1:5, 1:5)
#     camGrid = (1:5, 1:5)
#     θ = 0
#     flambda = 1
#     result = dualate(img, slmGrid, camGrid, θ, flambda)
#     @test size(result[1]) == (5, 5)

#     # Test 3: Check that the function works correctly with multiple images and a square grid
#     imgs = [rand(5, 5), rand(5, 5)]
#     result = dualate(imgs, slmGrid, camGrid, θ, flambda)
#     @test size(result[1]) == (5, 5)
#     @test size(result[2]) == (5, 5)

#     # # Test 4: Check that the function works correctly with a region of interest
#     # roi = (2:4, 2:4)
#     # result = dualate(img, slmGrid, camGrid, θ, flambda; roi=roi)
#     # @test size(result[1]) == (3, 3)
# end

