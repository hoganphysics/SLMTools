using Test
using SLMTools.LatticeTools.Resampling
using Interpolations: Flat


@testset "Resampling tests" begin
    @testset "downsample array" begin
        arr = reshape(1:16, (4, 4))
        @test isapprox(downsample(arr, 2), [3.5 11.5; 5.5 13.5], atol=1e-10)
    end

    @testset "downsample lattice" begin
        lattice = L = (1:12, 1:12)
        @test isapprox(downsample(lattice, 2)[1], 1.5:2.0:11.5, atol=1e-10)
        @test isapprox(downsample(lattice, 2)[2], 1.5:2.0:11.5, atol=1e-10)
    end

    @testset "upsample array" begin
        arr = collect(2:2:20)
        lat = (1:10,) .* 2
        newlat = (1:0.5:10,) .* 2
        uparr = upsample(arr, lat, newlat,bc=Flat()) .|> round
        @test uparr == 2.0:1.0:20.0
    end

end