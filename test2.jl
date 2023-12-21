using SLMTools
N= 1024
beam0 = lfGaussian(Intensity, (N, N), 5.0)

look(beam0.data)