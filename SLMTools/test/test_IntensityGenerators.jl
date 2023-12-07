using Test
using SLMTools.IntensityGenerators
using SLMTools.LatticeTools:natlat

@testset "IntensityGenerators tests" begin
    @testset "makeGaussian tests" begin
        sz = (10, 10)
        sigma = 1.0
        gaussian = makeGaussian(sz, sigma)

        @test size(gaussian) == sz
        @test isapprox(gaussian[1:3, 1:3], [0.082085 0.128735 0.182684; 0.128735 0.201897 0.286505; 0.182684 0.286505 0.40657], atol=1e-5)
    end

    @testset "makeRing tests" begin
        sz = tuple((rand(10:2:20) for _ in 1:2)...)
        r = 2
        w = 1.0
        ring = makeRing(sz, r, w)
        cind = sz .÷ 2
        center = getindex.(natlat(sz), cind)
        @test size(ring) == sz
        expected = exp(-(sqrt(center[1]^2 + center[2]^2) - r)^2 / (2 * w^2))
        @test isapprox(ring[cind...], expected)
    end

    @testset "makeLetterPattern tests" begin
        tw = 10
        NT = 100
        word = "AI"
        pattern = makeLetterPattern(tw, NT, word)

        @test size(pattern) == (NT, NT)
        # Add more tests based on the expected properties of the pattern
    end
end

