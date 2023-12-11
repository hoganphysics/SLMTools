using SLMTools, Plots, FFTW




n = 128
L0 = natlat((n, n))
f = LF{Generic}(rand(length.(L0)...), L0);


L1 = natlat((n,))
μ = LF{Intensity}(exp.(-L1[1] .^ 2 ./ 2), L1)
v = LF{Intensity}(exp.(-L1[1] .^ 2), L1);


as = 0.1:0.1:1.0
ϕs = [LF{RealPhase}(a * L1[1] .^ 2 / 2, L1) for a in as] |> Tuple;
divImgs = [sft(sqrt(μ) * ϕ) |> square for ϕ in ϕs] |> Tuple;

g = pdgs(divImgs, ϕs, 100, sqrt(v) * ϕs[1])

[μ.data, abs.(square(g).data)] |> plot

μ = LF{Intensity}([exp.(-(x^2 + y^2) ./ 2) for x in L0[1], y in L0[2]], L0)
v = LF{Intensity}([exp.(-(x^2 + y^2) ./ 4) for x in L0[1], y in L0[2]], L0)
p = gs(μ, v, 1000, LF{ComplexPhase}(zeros(size(μ)) * im, μ.L))
hcat([look(μ), look(sqrt(v)), look(p), look(sqrt(μ) * p)]...)