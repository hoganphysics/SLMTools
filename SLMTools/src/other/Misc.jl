module Misc
export ramp, nabs, centroid, window


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


end # module Misc