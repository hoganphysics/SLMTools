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
