module Misc
export ramp, nabs, centroid, window, normalizeDistribution, safeInverse, hyperSum


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
    centroid(img::Array{T,N}) where {T,N}

    Calculates the centroid of a multidimensional array `img`. The centroid is computed as a weighted average of indices, where the weights are the array elements.
    # Arguments
    - `img`: A multidimensional array.
    # Returns
    - A vector representing the centroid coordinates of `img`.
    """
centroid(img::Array{T,N}) where {T,N} = [sum((1:size(img)[j]) .* sum(img, dims=delete!(Set(1:N), j))[:]) for j = 1:N] ./ sum(img)

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







end # module Misc