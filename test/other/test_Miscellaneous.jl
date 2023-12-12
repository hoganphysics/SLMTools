using Test
using SLMTools

@testset "Miscellaneous tests" begin
    @testset "centroid" begin
        img1 = LF{Intensity}([1 2; 3 4], (1:2, 1:2))
        img2 = LF{Intensity}(ones(3,3), (1:3, 1:3))
        img3 = LF{Intensity}([1 0; 0 0], (1:2, 1:2))
        @test centroid(img1) == [1.7; 1.6]
        @test centroid(img1.data) == [1.7; 1.6]
        @test centroid(img2) == [2.0; 2.0]
        @test centroid(img3) == [1.0; 1.0]
    end

    @testset "window" begin
        # Test 1: Window of a 2D array
        img4 = ones(4, 4)
        w1 = 2
        @test size(window(img4, w1)) == (w1, w1)

        # Test 2: Window with odd window size
        w2 = 3
        @test size(window(img4, w2)) == (w2, w2)

        # Test 3: Window of a 3D array
        img5 = ones(4, 4, 4)
        w3 = 2
        @test size(window(img5, w3)) == (w3, w3, w3)
    end

    @testset "ramp tests" begin
        @test ramp(-1) == 0
        @test ramp(0) == 0
        @test ramp(1) == 1
    end

    @testset "nabs tests" begin
        v = [-1, 0, 1]
        @test nabs(v) == [1 / sqrt(2), 0, 1 / sqrt(2)]
    end

end

