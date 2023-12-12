module Misc
using ..LatticeTools
export ramp, nabs, centroid, window, normalizeDistribution, safeInverse, hyperSum, collapse, clip


"""
    ramp(x::T) where {T<:Number}

    Linear ramp function, returning zero for negative inputs and the input itself for non-negative inputs.

    # Arguments
    - `x::T`: A numerical input.

    # Returns
    - `0` if `x` is less than `0`, otherwise returns `x`.

    This function acts as a rectifier that clips all negative values to zero, which is useful in various mathematical and signal processing applications where negative values are not permissible.
    This function is also known in machine learning as the ReLU (Rectified Linear Unit) activation function.
    """
function ramp(x::T) where {T<:Number}
    return (x < 0 ? zero(T) : x)
end

"""
    nabs(v::AbstractArray{<:Number})

    Normalize the absolute values of a vector.

    # Arguments
    - `v::AbstractArray{<:Number}`: A vector of numerical values.

    # Returns
    - A vector where each element is the absolute value of the corresponding element in `v`, normalized by the Euclidean norm of `v`.

    This function can be used to normalize a vector in a context where the signs of its components are not relevant, only their magnitudes and relative proportions.
    """
nabs(v) = abs.(v) ./ sqrt(sum(abs.(v) .^ 2))



"""
    window(img::Array{T,N}, w) where {T,N}

    Creates a window around the centroid of an image `img`. The window has a size of `w` in each dimension.
    # Arguments
    - `img`: A multidimensional array representing an image.
    - `w`: The size of the window in each dimension.
    # Returns
    - Cartesian indices representing the window around the centroid of `img`.
    """
function window(img::Array{T,N}, w) where {T,N}
    c = CartesianIndex(round.(Int, centroid(img))...)
    return CartesianIndices(((1:w for i = 1:N)...,)) .- CartesianIndex((floor(Int, w / 2) for i = 1:N)...) .+ c
end

"""
    normalizeDistribution(U::AbstractArray{T,N}) where {T<:Number,N}

    Normalize an array to represent a probability distribution.

    This function takes an arbitrary numeric array `U` and transforms it into a positive real array 
    where the sum of all elements equals 1. This is commonly used to represent a discrete probability distribution.

    # Arguments
    - `U::AbstractArray{T,N}`: An array of any shape containing numeric values. The array can be of any dimension `N`.

    # Returns
    - `AbstractArray{T,N}`: An array of the same shape as `U`, where each element is non-negative, and 
    the sum of all elements is 1.
    """
function normalizeDistribution(U::AbstractArray{T,N}) where {T<:Number,N}
    # Makes an arbitrary array into a positive real array with sum 1, i.e. a probability distribution. 
    U = abs.(U)
    return U ./ sum(U)
end


"""
    safeInverse(x::Number)

    Calculate the multiplicative inverse of `x` safely.

    # Arguments
    - `x::Number`: The number to invert.

    # Returns
    - The multiplicative inverse of `x` if `x` is non-zero, otherwise returns zero of the same type as `x`.

    This function is a safeguard against division by zero errors by checking if `x` is zero before attempting to calculate its inverse.
    """
safeInverse(x::Number) = iszero(x) ? zero(typeof(x)) : 1 / x

"""
    hyperSum(A::AbstractArray{T,N}, originIdx::CartesianIndex, sumDim::Int, fixDims) where {T<:Number, N}

    Perform a cumulative sum of array `A` in a hyperplane fixed by `fixDims` at `originIdx`.

    # Arguments
    - `A::AbstractArray{T,N}`: An N-dimensional array to sum over.
    - `originIdx::CartesianIndex`: The index of the origin in the cumulative sum.
    - `sumDim::Int`: The dimension along which the cumulative sum is performed.
    - `fixDims`: The dimensions to be fixed (i.e., not summed over).

    # Returns
    - An array with the same dimensions as `A`, containing the cumulative sum over the specified dimension, adjusted by the value at `originIdx`.

    This function performs a cumulative sum over the dimension `sumDim` while fixing the other specified dimensions. It adjusts the sum relative to the value at `originIdx`, effectively setting the potential to zero at that point.
    """
function hyperSum(A::AbstractArray{T,N}, originIdx::CartesianIndex, sumDim::Int, fixDims) where {T<:Number,N}
    # Slices array A along dimensions fixDims and cumsums along dimension sumDim.  Does this for each value of 
    # the remaining dimensions.  originIdx indicates which CartesianIndex should be zero in the cumsum.
    sumDim in fixDims && error("sumDim can't be in fixDims")
    sumDim in 1:N && all(d in 1:N for d in fixDims) || error("More dimensions than array axes.")

    dims = (((i in fixDims || i == sumDim ? size(A)[i] : (originIdx[i]:originIdx[i])) for i = 1:N)...,)
    CI1 = CartesianIndices(dims)
    sumBox = cumsum(A[CI1], dims=sumDim)

    dims = (((i in fixDims ? size(A)[i] : (i == sumDim ? (originIdx[i]:originIdx[i]) : 1)) for i = 1:N)...,)
    CI2 = CartesianIndices(dims)
    originSlice = sumBox[CI2]
    return sumBox .- originSlice
end

"""
    collapse(x::AbstractArray{T,N}, i::Int) where {T<:Number, N}

    Sum all dimensions of an array except for the specified dimension.

    This function computes the sum of elements in the array `x` across all dimensions except the `i`-th dimension. The result is a collapse of the array into a lower-dimensional form based on the summation.

    # Arguments
    - `x::AbstractArray{T,N}`: The input array to be collapsed.
    - `i::Int`: The dimension that should not be summed.

    # Returns
    - An array where all dimensions except the `i`-th have been summed over.

    # Notes
    - The function preserves the `i`-th dimension while summing over the other dimensions.
    """
function collapse(x::AbstractArray{T,N}, i::Int) where {T<:Number,N}
    # Sums all dimensions of x except the i-th.  
    return sum(x, dims=(1:i-1..., i+1:N...))[:]
end

"""
    clip(x::Number, threshold::Number)

    Clip a number to zero if it falls below a specified threshold.

    This function evaluates a number `x` and sets it to zero if it is less than or equal to the specified `threshold`. If `x` is greater than the threshold, it is left unchanged.

    # Arguments
    - `x::Number`: The number to be evaluated.
    - `threshold::Number`: The threshold value for clipping.

    # Returns
    - The original number `x` if it is greater than `threshold`, otherwise `zero(x)`.

    # Notes
    - This function is useful for thresholding values, particularly in contexts where small numbers should be treated as zero for stability or simplicity.
    """
function clip(x::Number, threshold::Number)
    # Clips a number to zero if it's below the threshold. 
    x > threshold ? x : zero(x)
end

"""
    centroid(img::LF{Intensity,T,N}, threshold::Real=0.1) where {T, N}

    Calculate the centroid of a LatticeField in its own coordinates.

    This function computes the centroid (center of mass) of the intensity distribution in a LatticeField `img`. 
    Values below a specified `threshold` are clipped to zero before the calculation.

    # Arguments
    - `img::LF{Intensity,T,N}`: A LatticeField representing the intensity distribution.
    - `threshold::Real`: A threshold value for clipping intensities (default: 0.1).

    # Returns
    - A tuple representing the coordinates of the centroid within the LatticeField.

    # Errors
    - Throws an error if the total intensity sum is zero (i.e., a black image), as the centroid cannot be normalized in this case.

    The centroid is calculated by clipping the LatticeField intensities, summing them, and then computing 
    the weighted average position of the remaining intensity distribution.
    """
function centroid(img::LF{Intensity,T,N}, threshold::Real=0.1) where {T,N}
    # Centroid of a LatticeField, in its own coordinates. 
    clipped = clip.(img.data, threshold)
    s = sum(clipped)
    (s > 0) || error("Black image.  Can't normalize")
    return [(sum(collapse(clipped, i) .* img.L[i]) for i = 1:N)...] ./ s
end
centroid(img::Array{T,N}) where {T,N} = [sum((1:size(img)[j]) .* sum(img, dims=delete!(Set(1:N), j))[:]) for j = 1:N] ./ sum(img)



end # module Misc