module SLMTools


include("LatticeTools/LatticeTools.jl")
using .LatticeTools
export Lattice, elq, RealPhase, Generic, Phase, FieldVal, ComplexPhase, UPhase, UnwrappedPhase, S1Phase, Intensity, Amplitude, Modulus, RealAmplitude, RealAmp, ComplexAmplitude, ComplexAmp, LatticeField, LF, subfield, wrap, square, sublattice, normalizeLF, phasor
export natrange, natlat, sft, isft, padout, naturalize, latticeDisplacement, toDim, r2, ldot, Nyquist
export downsample, upsample, coarsen, sublattice
export dualLattice, dualShiftLattice, ldq, dualPhase

include("LFIO/Misc.jl")
using .Misc
export ramp, nabs, window, normalizeDistribution, safeInverse, hyperSum, centroid, clip, collapse, SchroffError

include("LFIO/ImageProcessing.jl")
using .ImageProcessing
export getImagesAndFilenames, imageToFloatArray, itfa, castImage, loadDir, parseFileName, parseStringToNum, getOrientation, dualate, linearFit, savePhase, saveBeam
# saveAs8BitBMP

include("LFIO/LFTemplates.jl")
using .LFTemplates
export lfParabola, lfGaussian, lfRing, lfCap, ftaText, lfText, lfRect, lfRand

include("PhaseRetrieval/IFT.jl")
using .IFT
export gs, gsLog, gsIter, pdgs, pdgsIter, pdgsLog, oneShot, pdgsError, mraf

include("PhaseRetrieval/OT.jl")
using .OT
export getCostMatrix, pdCostMatrix, pdotBeamEstimate,  mapify, scalarPotentialN, otPhase, pdotPhase, pdotBeamEstimate

include("LFIO/Visualization.jl")
using .Visualization
export look


export testing
end # module SLMTools
