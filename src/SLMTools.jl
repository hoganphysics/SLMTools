module SLMTools

using Interpolations: Flat, Periodic, CubicSplineInterpolation, LinearInterpolation, Linear
using FFTW  # Needed for IFT.jl
#using FFTW: fftshift, ifftshift, fft, ifft, plan_fft, plan_ifft  # Needed for DualLattices.jl and IFT.jl
using LinearAlgebra 
using ToeplitzMatrices 
using OptimalTransport: sinkhorn
#using LinearAlgebra: norm
using FileIO: load, save
using Images: Gray, RGB, Colorant  # Needed for ImageProcessing.jl
#using Images: Gray  # Needed for Visualization.jl
using FreeTypeAbstraction: findfont, renderstring!
using ColorTypes
using FixedPointNumbers   # Gray, N0f8

include("LatticeFields/LatticeField.jl")

include("LatticeFields/LatticeUtils.jl")

include("LatticeFields/Resampling.jl")

include("LatticeFields/DualLattices.jl")

include("LFIO/Misc.jl")

include("LFIO/ImageProcessing.jl")

include("LFIO/LFTemplates.jl")

include("PhaseRetrieval/IFT.jl")

include("PhaseRetrieval/OT.jl")

include("LFIO/Visualization.jl")

include("LFIO/BMP8Writer.jl")

end 
