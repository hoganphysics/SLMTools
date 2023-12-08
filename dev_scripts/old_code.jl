"""
    getMoment(A::AbstractArray{<:Number,N}, m::NTuple{M,Integer}) where {N,M}

    Compute the moment of a given array `A`.

    # Arguments
    - `A::AbstractArray{<:Number,N}`: A multidimensional array representing a discrete function or data set.
    - `m::NTuple{M,Integer}`: A tuple specifying the indices for the moment calculation. For instance, `(i_1, i_2, ..., i_M)` indicates that the moment `<x_{i_1} x_{i_2} ... x_{i_M}>` should be calculated.

    # Returns
    - The calculated moment of the array `A` based on the specified indices in `m`.

    # Details
    This function calculates the moment of the array `A` along the specified dimensions in `m`. It internally handles non-dimensional and dimensional aspects of the array to focus the moment calculation on the relevant axes.
    """
function getMoment(A::AbstractArray{<:Number,N}, m::NTuple{M,Integer}) where {N,M}
    # Returns the moment <x_{i_1} x_{i_2} ... x_{i_M}> of A, where m = (i_1,i_2,...,i_M)
    pows = zeros(Int, N)
    for j in m
        pows[j] += 1
    end
    nondims = filter(x -> pows[x] == 0, 1:N)
    indims = filter(x -> pows[x] != 0, 1:N)
    B = sum(A, dims=nondims)
    return sum(B .* [prod([I[j]^pows[j] for j in indims]) for I in CartesianIndices(size(B))])
end


#TODO do i need this?
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

#TODO do i need this?
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