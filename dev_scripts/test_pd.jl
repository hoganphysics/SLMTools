using SLMTools
using Test

n = 128
L0 = natlat((n, n))
L1 = natlat((n,))
μ = LF{Intensity}(exp.(-L1[1] .^ 2 ./ 2), L1)
v = LF{Intensity}(exp.(-L1[1] .^ 2), L1);


as = 0.1:0.1:1.0
ϕs = [LF{RealPhase}(a * L1[1] .^ 2 / 2, L1) for a in as] |> Tuple;
divImgs = [sft(sqrt(μ) * ϕ) |> square for ϕ in ϕs] |> Tuple;

g = pdgs(divImgs, ϕs, 100, sqrt(v) * ϕs[1])
μ = LF{Intensity}([exp.(-(x^2 + y^2) ./ 2) for x in L0[1], y in L0[2]], L0)
v = LF{Intensity}([exp.(-(x^2 + y^2) ./ 4) for x in L0[1], y in L0[2]], L0)
p = gs(μ, v, 1000, LF{ComplexPhase}(zeros(size(μ)) * im, μ.L))

@testset "phase diversity tests" begin
    @testset "gs algorithm alone" begin
        @test isa(g, LatticeField{ComplexAmplitude})
        @test g[48:50].data ≈ [0.5510345631134064 + 0.14050740448745414im
            0.5877216316040312 + 0.14985954977667548im
            0.6244073206579513 + 0.1592112010176915im]
    end    
    @testset "pdgs algorithm alo" begin
        @test isapprox(p[48:50, 48:50].data, [0.143712-0.989619im -0.0419511-0.99912im -0.727503-0.686104im
                0.919021-0.394209im 0.659497-0.751708im -0.0029928-0.999996im
                0.957521-0.288362im 0.76654-0.642196im 0.13438-0.99093im], atol=1e-5)
        @test isa(p, LatticeField{ComplexPhase})
        @test all(isapprox.(p.L, (-5.65685424949238:0.08838834764831843:5.568465901844061, -5.65685424949238:0.08838834764831843:5.568465901844061)))
    end
    
end

# hcat([look(μ), look(sqrt(v)), look(p), look(sqrt(μ) * p)]...)