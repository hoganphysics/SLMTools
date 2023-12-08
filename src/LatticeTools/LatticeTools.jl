module LatticeTools

include("LatticeCore.jl")
using .LatticeCore
export Lattice, natlat, sft, isft, padout

include("Resampling.jl")
using .Resampling
export downsample, upsample, coarsen, sublattice

include("DualLattices.jl")
using .DualLattices
export dualShiftLattice, dualLattice, dualate


# export  sublattice, latticeDisplacement
# export toDim, shiftedPhases, dualPhase



"""
    latticeDisplacement(L::Lattice)

    Calculates the displacement of each dimension in a lattice `L`. The displacement is computed as the midpoint offset of each lattice dimension.

    # Arguments
    - `L`: A lattice.

    # Returns
    - A vector of displacements for each dimension of the lattice.
    """
function latticeDisplacement(L::Lattice)
    return [l[1] + floor(length(l) / 2) * step(l) for l in L]
end

"""
    toDim(v, d::Int, n::Int)

    Reshapes a vector `v` into an `n`-dimensional array where all dimensions except the `d`-th are of size 1.

    # Arguments
    - `v`: A vector to be reshaped.
    - `d`: The dimension to retain the length of `v`.
    - `n`: The number of dimensions of the resulting array.

    # Returns
    - An `n`-dimensional array with the `d`-th dimension equal to the length of `v` and other dimensions of size 1.
    """
function toDim(v, d::Int, n::Int)
    reshape(v, ((i == d ? length(v) : 1) for i = 1:n)...)
end

"""
    dualPhase(L::Lattice{N}) where {N}

    Computes the dual phase of a lattice `L`. This function involves shifting the lattice, computing its displacement, and applying a complex exponential transformation.

    # Arguments
    - `L`: A lattice.

    # Returns
    - An array representing the dual phase of the lattice `L`.
    """
function dualPhase(L::Lattice{N}) where {N}
    dL = dualShiftLattice(L, 1)
    b = latticeDisplacement(L)
    return exp.(-2 * pi * im .* .+((b[i] * toDim(dL[i], i, N) for i = 1:N)...))
end

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
