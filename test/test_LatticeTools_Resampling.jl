using Test
using SLMTools.LatticeTools.Resampling

@testset "Resampling tests" begin
    @testset "downsample array" begin
        arr = reshape(1:16, (4, 4))
        @test downsample(arr, 2) == [14 46; 22 54]
    end

    @testset "downsample lattice" begin
        lattice = L = (1:12, 1:12)
        @test downsample(lattice, 2) == (1.5:2:12.5, 1.5:2:12.5)
    end

    @testset "upsample array" begin
        arr = collect(2:2:20)
        lat = (1:10,) .* 2
        newlat = (1:0.5:10,) .* 2
        uparr = upsample(arr, lat, newlat) .|> round
        @test uparr == 2.0:1.0:20.0
    end

    @testset "upsampleBeam" begin
        beam = [1+1im 2+2im; 3+3im 4+4im] .* 1.0
        @test upsampleBeam(beam, 2) ≈ [4.0+4.0im 3.5+3.5im 4.0+4.0im 3.5+3.5im
            3.0+3.0im 2.5+2.5im 3.0+3.0im 2.5+2.5im
            4.0+4.0im 3.5+3.5im 4.0+4.0im 3.5+3.5im
            3.0+3.0im 2.5+2.5im 3.0+3.0im 2.5+2.5im]
    end
end