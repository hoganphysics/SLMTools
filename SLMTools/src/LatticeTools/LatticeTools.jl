module LatticeTools

include("LatticeCore.jl")
using .LatticeCore: Lattice, natlat
export Lattice, natlat, sft, isft

include("Resampling.jl")
using .Resampling
export downsample, upsample, upsampleBeam

include("DualLattices.jl")
using .DualLattices
export dualShiftLattice, dualLattice, dualate


# export  sublattice, latticeDisplacement
# export toDim, shiftedPhases, dualPhase





"""
    sublattice(L::Lattice{N}, box::CartesianIndices) where {N}

Extracts a sublattice from a given `Lattice` `L` using the specified `box` of Cartesian indices.

# Arguments
- `L`: A lattice from which the sublattice is to be extracted.
- `box`: Cartesian indices defining the region of the sublattice.

# Returns
- A tuple representing the extracted sublattice.
- Throws an error if the dimensions of the `box` do not match the lattice.
"""
function sublattice(L::Lattice{N}, box::CartesianIndices) where {N}
    length(size(box)) == N || error("box should have same number of dimensions as lattice.")
    return ((L[j][box[1][j]:box[end][j]] for j = 1:N)...,)
end

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


# """
#     downsample(μ::AbstractArray{T,N}, r::Int; reducer=sum) where {T, N}

#     Downsampling operation on an N-dimensional array `μ` by a factor of `r` using a specified reduction function (default is sum). This function reduces the size of the array by grouping and aggregating elements.

#     # Arguments
#     - `μ`: An N-dimensional array to be downsampled.
#     - `r`: Downsampling factor, which must evenly divide each dimension of `μ`.
#     - `reducer`: A function to aggregate elements in each group (default is `sum`).

#     # Returns
#     - A downsampled version of the input array `μ`.

#     # Errors
#     - Throws an error if `r` does not evenly divide each dimension of `μ`.
# """
# function downsample(μ::AbstractArray{T,N}, r::Int; reducer=sum) where {T,N}
#     size(μ) .% r == 0 .* size(μ) || error("Downsample factor must divide array size.")
#     CR = CartesianIndices(((0:r-1 for i = 1:N)...,))
#     CD = CartesianIndices(((1:r:s for s in size(μ))...,))
#     return [reducer(μ[I.+CR]) for I in CD]
# end

# """
#     downsample(L::Lattice{N}, r::Int) where {N}

#     Downsampling operation for a lattice `L` by a factor of `r`. This function reduces the number of points in the lattice by selecting points at regular intervals based on the downsampling factor.

#     # Arguments
#     - `L`: A lattice to be downsampled.
#     - `r`: Downsampling factor, which must evenly divide the length of each dimension of `L`.

#     # Returns
#     - A downsampled version of the lattice `L`.

#     # Errors
#     - Throws an error if `r` does not evenly divide the length of each dimension of `L`.
# """
# function downsample(L::Lattice{N}, r::Int) where {N}
#     length.(L) .% r == 0 .* length.(L) || error("Downsample factor must divide lattice size.")
#     return ((((0:r:length(l)-1) .+ (r - 1) / 2) .* step(l) .+ l[1] for l in L)...,)
# end

# """
#     upsample(A::AbstractArray{T,N}, LL::Lattice{N}, LH::Lattice{N}) where {T, N}

#     Upsamples an N-dimensional array `A` from a lower-resolution lattice `LL` to a higher-resolution lattice `LH` using cubic spline interpolation.

#     # Arguments
#     - `A`: An N-dimensional array to be upsampled.
#     - `LL`: The lattice corresponding to the lower-resolution of `A`.
#     - `LH`: The lattice defining the higher-resolution target.

#     # Returns
#     - An upsampled version of the array `A` fitting the higher-resolution lattice `LH`.

#     # Note
#     - This function uses cubic spline interpolation for the upsampling process.
# """
# function upsample(A::AbstractArray{T,N}, LL::Lattice{N}, LH::Lattice{N}) where {T,N}
#     return cubic_spline_interpolation(LL, A, extrapolation_bc=Flat())[LH...]
# end

# """
#     upsampleBeam(beam::Matrix{ComplexF64}, st::Int)

#     Upsamples a complex beam matrix using cubic spline interpolation. This is typically used to increase the resolution of a beam's spatial representation.

#     # Arguments
#     - `beam`: A matrix representing the complex field of a beam.
#     - `st`: Step size for upsampling.

#     # Returns
#     - A matrix representing the upsampled beam.
# """
# function upsampleBeam(beam::Matrix{ComplexF64}, st::Int)
#     grid = (((1:st:st*s) .- 1 for s in size(beam))...,)
#     return cubic_spline_interpolation(grid, ifftshift(beam), extrapolation_bc=Periodic())[CartesianIndices(((0:st*s-1 for s in size(beam))...,))] |> fftshift


# end









end # module LatticeTools
