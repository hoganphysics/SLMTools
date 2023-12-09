module LatticeTools

include("LatticeCore.jl") # uses LatticeFields.jl
using .LatticeCore
export Lattice, elq, RealPhase, Generic, Phase, FieldVal, ComplexPhase, UPhase, UnwrappedPhase, S1Phase, Intensity, Amplitude, Modulus, RealAmplitude, RealAmp, ComplexAmplitude, ComplexAmp, LatticeField, LF, subfield, wrap, square, sublattice
export natlat, sft, isft, padout

include("Resampling.jl")
using .Resampling
export downsample, upsample, coarsen, sublattice

include("DualLattices.jl")
using .DualLattices
export dualShiftLattice, dualLattice, dualate, dualPhase, ldq


export   latticeDisplacement, toDim







"""
    shiftedPhases(lattices::Vector{S}, phases::Vector{Array{R,N}}) where {S<:Lattice{N}, R<:Complex, N}

    Applies a dual phase shift to a set of phases corresponding to a vector of lattices. This is used to modify the phase of light in each lattice in a complex way.

    # Arguments
    - `lattices`: A vector of lattice structures.
    - `phases`: A vector of complex arrays representing the phases corresponding to each lattice.

    # Returns
    - A vector of complex arrays, each representing the phase-shifted version of the corresponding input phase.
    """
function shiftedPhases(lattices::Vector{S}, phases::Vector{Array{R,N}}) where {S<:Lattice{N},R<:Complex} where {N}
    return [phases[j] .* dualPhase(lattices[j]) for j = 1:length(lattices)]
end





end # module LatticeTools
