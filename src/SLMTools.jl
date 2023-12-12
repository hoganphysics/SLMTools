module SLMTools


include("LatticeTools/LatticeTools.jl")
using .LatticeTools
export Lattice, elq, RealPhase, Generic, Phase, FieldVal, ComplexPhase, UPhase, UnwrappedPhase, S1Phase, Intensity, Amplitude, Modulus, RealAmplitude, RealAmp, ComplexAmplitude, ComplexAmp, LatticeField, LF, subfield, wrap, square, sublattice
export natlat, sft, isft, padout, naturalize, latticeDisplacement, toDim
export downsample, upsample, coarsen, sublattice
export dualLattice, dualShiftLattice, dualate, ldq, dualPhase


include("other/ImageProcessing.jl")
using .ImageProcessing
export getImagesAndFilenames, imageToFloatArray, itfa, getCamGrid, sweepCoords, findCOM, lineFit, castImage, loadDir, parseFileName, parseStringToNum, getOrientation

include("other/PhaseIntensityMasks.jl")
using .PhaseIntensityMasks
export lfRampedParabola, lfGaussian, lfRing

include("other/Misc.jl")
using .Misc
export ramp, nabs, centroid, window, normalizeDistribution, safeInverse, hyperSum

include("other/HoganParameters.jl")
using .HoganParameters
export dxcam, dxslm, nslm, flambda, dXslm, Lslm, dL, LslmE, dLE

include("PhaseDiversity/LatticeIFTA.jl")
using .LatticeIFTA
export phasor, gs, gsIter, pdgs, pdgsIter, oneShot

include("PhaseDiversity/OTHelpers.jl")
export getCostMatrix, pdCostMatrix, pdotBeamEstimate,  mapify, scalarPotentialN,  OTphase, pdotPhase

include("other/VisualizationHelpers.jl")
using .VisualizationHelpers
export look

include("other/FileHelpers.jl")
export savePhase, saveBeam, saveAs8BitBMP


end # module SLMTools