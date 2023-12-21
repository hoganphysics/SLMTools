using SLMTools
# saveAs8BitBMP(convert.(Int64, floor.(rand(100,100).*255)), "test.bmp")
using Plots
N = 1024
NSLM = 1024
sz = (10, 10)
sigma = 0.5
beam0 = lfGaussian(Intensity, (N, N), 5.0)
dsbeam0 = downsample(beam0, 8) # Make a 8x downsampled image
p1 = look(dsbeam0)

target = lfRing(Intensity, (N, N), 5, 1)
p2 = look(target)
dstarget = downsample(target, 8)
p3 = look(dstarget)

L = natlat((128, 128))
phi1 = otPhase(dsbeam0, dstarget, 0.1)
p4 = look(phi1)

phi2 = upsample(phi1, 8)
p5 = look(phi2)

phi = phi2
beam_out = sft(sqrt.(dsbeam0.data) .* exp.(2pi * im * phi1.data)) .^ 2 |> nabs;
# p6 = look(beam_out)
LFtest = LatticeField{Intensity}(beam_out, dsbeam0.L);



# phi = phi2.data
# # phi.+= 2 * pi
# phiout = convert(Matrix{Int64}, floor.(((phi .+ abs(minimum(phi))) .% (2 * pi) .+ 0.0) .* 256 / (2 * pi)))
# minimum(phiout), maximum(phiout)

# # saveAs8BitBMP(phiout, "ringplan.bmp")

# # plot(1:10,1:10)