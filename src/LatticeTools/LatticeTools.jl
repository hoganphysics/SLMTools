module LatticeTools

include("LatticeCore.jl") # uses LatticeFields.jl
using .LatticeCore
export Lattice, elq, RealPhase, Generic, Phase, FieldVal, ComplexPhase, UPhase, UnwrappedPhase, S1Phase, Intensity, Amplitude, Modulus, RealAmplitude, RealAmp, ComplexAmplitude, ComplexAmp, LatticeField, LF, subfield, wrap, square, sublattice
export natlat, sft, isft, padout, latticeDisplacement, toDim

include("Resampling.jl")
using .Resampling
export downsample, upsample, coarsen, sublattice

include("DualLattices.jl")
using .DualLattices
export dualShiftLattice, dualLattice, dualate, dualPhase, ldq




end # module LatticeTools
