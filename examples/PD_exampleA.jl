using SLMTools, Plots, FFTW

sft(v) = fftshift(fft(ifftshift(v)))
isft(v) = fftshift(ifft(ifftshift(v)))
sft(v::LF{ComplexAmplitude}) = LF{ComplexAmplitude}(fftshift(fft(ifftshift(v.data))), dualShiftLattice(v.L, v.flambda), v.flambda)
isft(v::LF{ComplexAmplitude}) = LF{ComplexAmplitude}(fftshift(ifft(ifftshift(v.data))), dualShiftLattice(v.L, v.flambda), v.flambda)


n = 128
L0 = natlat((n, n))
f = LF{Generic}(rand(length.(L0)...), L0);


L1 = natlat((n,))
μ = LF{Intensity}(exp.(-L1[1] .^ 2 ./ 2), L1)
ν = LF{Intensity}(exp.(-L1[1] .^ 2), L1);


αs = 0.1:0.1:1.0
ϕs = [LF{RealPhase}(α * L1[1] .^ 2 / 2, L1) for α in αs] |> Tuple;
divImgs = [sft(sqrt(μ) * ϕ) |> square for ϕ in ϕs] |> Tuple;

g = pdgs(divImgs, ϕs, 100, sqrt(ν) * ϕs[1])

[μ.data, abs.(square(g).data)] |> plot

μ = LF{Intensity}([exp.(-(x^2 + y^2) ./ 2) for x in L0[1], y in L0[2]], L0)
ν = LF{Intensity}([exp.(-(x^2 + y^2) ./ 4) for x in L0[1], y in L0[2]], L0)
p = gs(μ, ν, 1000, LF{ComplexPhase}(zeros(size(μ)) * im, μ.L))
# look(p)
# hcat([look(μ), look(sqrt(ν)), look(p), look(sqrt(μ) * p)]...)