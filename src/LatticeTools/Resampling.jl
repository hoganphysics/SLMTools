module Resampling
using ..LatticeCore: Lattice
using Interpolations: cubic_spline_interpolation, Flat, Periodic
using FFTW: fftshift, ifftshift

export downsample, upsample, upsampleBeam, sublattice
export upsample, upsampleBeam

#region downsample
"""
    downsample(μ::AbstractArray{T,N}, r::Int; reducer=sum) where {T, N}

    Downsampling operation on an N-dimensional array `μ` by a factor of `r` using a specified reduction function (default is sum). This function reduces the size of the array by grouping and aggregating elements.

    # Arguments
    - `μ`: An N-dimensional array to be downsampled.
    - `r`: Downsampling factor, which must evenly divide each dimension of `μ`.
    - `reducer`: A function to aggregate elements in each group (default is `sum`).

    # Returns
    - A downsampled version of the input array `μ`.

    # Errors
    - Throws an error if `r` does not evenly divide each dimension of `μ`.
"""
function downsample(μ::AbstractArray{T,N}, r::Int; reducer=sum) where {T,N}
    size(μ) .% r == 0 .* size(μ) || error("Downsample factor must divide array size.")
    CR = CartesianIndices(((0:r-1 for i = 1:N)...,))
    CD = CartesianIndices(((1:r:s for s in size(μ))...,))
    return [reducer(μ[I.+CR]) for I in CD]
end

"""
    downsample(L::Lattice{N}, r::Int) where {N}

    Downsampling operation for a lattice `L` by a factor of `r`. This function reduces the number of points in the lattice by selecting points at regular intervals based on the downsampling factor.

    # Arguments
    - `L`: A lattice to be downsampled.
    - `r`: Downsampling factor, which must evenly divide the length of each dimension of `L`.

    # Returns
    - A downsampled version of the lattice `L`.

    # Errors
    - Throws an error if `r` does not evenly divide the length of each dimension of `L`.
"""
function downsample(L::Lattice{N}, r::Int) where {N}
    length.(L) .% r == 0 .* length.(L) || error("Downsample factor must divide lattice size.")
    return ((((0:r:length(l)-1) .+ (r - 1) / 2) .* step(l) .+ l[1] for l in L)...,)
end
#endregion

#region upsample
"""
    upsample(A::AbstractArray{T,N}, LL::Lattice{N}, LH::Lattice{N}) where {T, N}

    Upsamples an N-dimensional array `A` from a lower-resolution lattice `LL` to a higher-resolution lattice `LH` using cubic spline interpolation.

    # Arguments
    - `A`: An N-dimensional array to be upsampled.
    - `LL`: The lattice corresponding to the lower-resolution of `A`.
    - `LH`: The lattice defining the higher-resolution target.

    # Returns
    - An upsampled version of the array `A` fitting the higher-resolution lattice `LH`.

    # Note
    - This function uses cubic spline interpolation for the upsampling process.
"""
function upsample(A::AbstractArray{T,N}, LL::Lattice{N}, LH::Lattice{N}) where {T,N}
    return cubic_spline_interpolation(LL, A, extrapolation_bc=Flat())[LH...]
end

"""
    upsampleBeam(beam::Matrix{ComplexF64}, st::Int)

    Upsamples a complex beam matrix using cubic spline interpolation. This is typically used to increase the resolution of a beam's spatial representation.

    # Arguments
    - `beam`: A matrix representing the complex field of a beam.
    - `st`: Step size for upsampling.

    # Returns
    - A matrix representing the upsampled beam.
"""
function upsampleBeam(beam::Matrix{ComplexF64}, st::Int)
    grid = (((1:st:st*s) .- 1 for s in size(beam))...,)
    return cubic_spline_interpolation(grid, ifftshift(beam), extrapolation_bc=Periodic())[CartesianIndices(((0:st*s-1 for s in size(beam))...,))] |> fftshift
end


"""
    upsample(f::LatticeField{S,T,N}, Lu::Lattice{N}; interpolation=cubic_spline_interpolation, bc=Periodic()) where {S<:FieldVal,T<:Number,N}

Upsamples a `LatticeField` instance to a new lattice `Lu`. This function interpolates the data from the 
original lattice of `f` to the new lattice `Lu` using the specified interpolation method and boundary conditions.

# Parameters
- `f`: An instance of `LatticeField` to be upsampled.
- `Lu`: The new lattice to which the field will be upsampled.
- `interpolation`: The interpolation method used for upsampling. Defaults to `cubic_spline_interpolation`.
- `bc`: Boundary conditions for the interpolation. Defaults to `Periodic()`.

# Returns
A new `LatticeField` instance with data upsampled to the lattice `Lu`.

The `upsample` function is versatile, allowing for the interpolation of `f` onto any lattice `Lu`, regardless of the 
relationship between the original and new lattices. It is particularly useful in applications where adjusting the 
resolution or grid spacing of data is required, such as in image processing or numerical simulations in optical trapping.
"""
function upsample(x::AbstractArray{T,N}, Ld::Lattice{N}, Lu::Lattice{N}; interpolation=cubic_spline_interpolation, bc=Periodic()) where {T<:Number,N}
    # Upsamples an array from coarse lattice Ld to fine lattice Lu.  Actually, there doesn't need to be any relationship
    # between the lattices at all--this function just interpolates x from one lattice to the other. 
    if size(x) != length.(Ld)
        throw(DomainError((size(x), length.(Ld)), "upsample: Size of array x does not match size of lattice Ld."))
    end
    return interpolation(Ld, x, extrapolation_bc=bc)[Lu...]
end


#endregion

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

function sublattice(L::Lattice{N}, x::Union{I,AbstractRange{I}}...) where {N,I<:Integer}
    length(L) != length(x) && throw(DimensionMismatch("Wrong number of indices while attempting to make sublattice."))
    y = (((i isa Integer) ? (i:i) : i for i in x)...,)
    return ((L[j][y[j]] for j = 1:N)...,)
end

end