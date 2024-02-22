using SLMTools
N= 1024
saveAs8BitBMP(convert.(Int64, floor.(rand(100, 100) .* 255)), "test.bmp")
beam0 = lfGaussian(Intensity, (N, N), 5.0)

look(beam0.data)