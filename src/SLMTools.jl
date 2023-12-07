module SLMTools
# using FileIO, Plots, Images, FreeTypeAbstraction
# using Interpolations, FFTW, Statistics, OptimalTransport
# using PyCall

# include("SubModule/Submodule.jl")
# using .Submodule
# export subtestA, subtestB, subsubtestA, subsubtestB

include("LatticeTools/LatticeTools.jl")
using .LatticeTools
export Lattice, natlat, sft, isft
export downsample, upsample, upsampleBeam
export dualLattice, dualShiftLattice, dualate
# export toDim, shiftedPhases, dualPhase,dualShiftLattice, padout, dualate, dualLattice, sft, isft, centroid, window, sublattice, latticeDisplacement

include("other/ImageProcessing.jl")
using .ImageProcessing
export getImagesAndFilenames, imageToFloatArray, itfa, getCamGrid, sweepCoords, findCOM, lineFit

include("masks/PhaseGenerators.jl")
using .PhaseGenerators
export makeRampedParabola

include("masks/IntensityGenerators.jl")
using .IntensityGenerators
export makeGaussian, makeRing, makeLetter, makeLetterPattern

include("other/Misc.jl")
using .Misc
export ramp, nabs, centroid, window

include("other/HoganParameters.jl")
using .HoganParameters
export dxcam, dxslm, nslm, flambda, dXslm, Lslm, dL, LslmE, dLE

# include("EstimateInverse.jl")
# export estimateInverseField, estimateInverseImage

# include("GSalgorithms.jl")
# export gs, pdgs, ampSwap, phasor, pdError1, cpdError1, croppedPDGS, arraySizeChecks, gsRefine

# include("OTHelpers.jl")
# export getCostMatrix, pdCostMatrix, pdotBeamEstimate, safeInverse, mapify, scalarPotential2, scalarPotentialN, hyperSum, OTphase, OTreducedPhase

# include("VisualizationHelpers.jl")
# export plotBeamMap, hm

# include("FileHelpers.jl")
# export savePhase, saveBeam, saveAs8BitBMP


end # module SLMTools