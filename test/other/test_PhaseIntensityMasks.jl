using Test
# using SLMTools.PhaseGenerators
using SLMTools


@testset "mask tests" begin
    arr0 = lfParabola(RealPhase, natlat(10, 10), 5 / sqrt(2), (5, 5), flambda=1, center=(1, 0))
    arr = arr0.data
    @testset "RampedParabola tests" begin
        @test arr[5, :] ≈ [6.40000000e-01, 3.60000000e-01, 1.60000000e-01, 4.00000000e-02, 0.00000000e+00, 4.00000000e-02, 1.60000000e-01, 3.60000000e-01, 6.40000000e-01, 2.22044605e-16]
        @test arr[:, 5] ≈ [0.5086291501015239, 0.511471862576143, 0.5943145750507619, 0.7571572875253809, 0.0, 0.32284271247461904, 0.7256854249492382, 0.20852813742385723, 0.7713708498984764, 0.41421356237309537]
    end
end
@testset "lfGaussian Tests" begin
    sz = natlat(10, 10)
    sigma = 0.5
    lf = lfGaussian(Intensity, sz, sigma)

    @test typeof(lf) == LF{Intensity,Float64,2}
    @test size(lf.data) == sz
    @test lf.flambda == 1.0
    @test all(lf.data .<= 1.0)
end
@testset "lfRing Tests" begin
    sz = natlat(10, 10)
    r = 3.0
    w = 1.0
    lf = lfRing(Intensity, sz, r, w)

    @test typeof(lf) == LF{Intensity,Float64,2}
    @test size(lf.data) == size(sz)
    @test lf.flambda == 1.0
    @test all(lf.data .<= 1.0)
end

# @testset "makeLetterPattern tests" begin
#     tw = 10
#     NT = 100
#     word = "AI"
#     pattern = makeLetterPattern(tw, NT, word)

#     @test size(pattern) == (NT, NT)
#     # Add more tests based on the expected properties of the pattern
# end