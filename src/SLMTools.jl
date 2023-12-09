module SLMTools


include("LatticeTools/LatticeTools.jl")
using .LatticeTools
export Lattice, elq, RealPhase, Generic, Phase, FieldVal, ComplexPhase, UPhase, UnwrappedPhase, S1Phase, Intensity, Amplitude, Modulus, RealAmplitude, RealAmp, ComplexAmplitude, ComplexAmp, LatticeField, LF, subfield, wrap, square, sublattice
export natlat, sft, isft, padout
export downsample, upsample, coarsen, sublattice
export dualLattice, dualShiftLattice, dualate, ldq, dualPhase
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

include("PhaseDiversity/LatticeIFTA.jl")
using .LatticeIFTA
export phasor, gs, gsIter, pdgs, pdgsIter, oneShot

# include("OTHelpers.jl")
# export getCostMatrix, pdCostMatrix, pdotBeamEstimate, safeInverse, mapify, scalarPotential2, scalarPotentialN, hyperSum, OTphase, OTreducedPhase

# include("VisualizationHelpers.jl")
# export plotBeamMap, hm

# include("FileHelpers.jl")
# export savePhase, saveBeam, saveAs8BitBMP


end # module SLMTools