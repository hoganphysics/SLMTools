module PhaseIntensityMasks

using ..LatticeTools

export lfRampedParabola, lfGaussian, lfRing

"""
    lfRampedParabola(T::DataType, imSize::NTuple{N,Integer}, pixelScale::Real, center::NTuple{N,Real}, 
                     pAmp::Real, rAmps::NTuple{N,Real}; L::Union{Nothing, Lattice{N}}=nothing, flambda::Real=1) where N

    Create a LatticeField representing a ramped parabolic distribution.

    This function generates a parabolic intensity distribution with additional linear ramps. The distribution is 
    defined on a lattice and wrapped to be within [0,1).

    # Arguments
    - `T::DataType`: The data type of the LatticeField.
    - `imSize::NTuple{N,Integer}`: The size of the image to be generated.
    - `pixelScale::Real`: The scaling factor for the pixel size.
    - `center::NTuple{N,Real}`: The center of the parabola.
    - `pAmp::Real`: The amplitude of the parabola.
    - `rAmps::NTuple{N,Real}`: The amplitudes of the linear ramps in each dimension.
    - `L::Union{Nothing, Lattice{N}}`: The lattice structure (default: natural lattice).
    - `flambda::Real`: Scaling factor for the wavelength (default: 1).

    # Returns
    - `LF{T}`: A LatticeField with the ramped parabolic distribution.
    """
function lfRampedParabola(T::DataType, imSize::NTuple{N,Integer}, pixelScale::Real, center::NTuple{N,Real},     pAmp::Real, rAmps::NTuple{N,Real}; L::Union{Nothing,Lattice{N}}=nothing, flambda::Real=1) where {N}
    if isnothing(L)
        L = natlat(imSize)
    end
    im = ([pAmp * sum((Tuple(I) .- center) .^ 2) / (2 * pixelScale^2) for I in CartesianIndices(imSize)] .+
          [sum((Tuple(I) .- center) .* rAmps) / pixelScale for I in CartesianIndices(imSize)]) |> (x -> mod.(x, 1))
    return LF{T}(im, L, flambda)
end

"""
    lfGaussian(T::DataType, sz::NTuple{N,Integer}, sigma::Number; L::Union{Nothing, Lattice{N}}=nothing, flambda::Real=1) where N

    Create a LatticeField representing a Gaussian distribution.

    This function generates a Gaussian intensity distribution, defined on a lattice.

    # Arguments
    - `T::DataType`: The data type of the LatticeField.
    - `sz::NTuple{N,Integer}`: The size of the Gaussian distribution.
    - `sigma::Number`: The standard deviation of the Gaussian distribution.
    - `L::Union{Nothing, Lattice{N}}`: The lattice structure (default: natural lattice).
    - `flambda::Real`: Scaling factor for the wavelength (default: 1).

    # Returns
    - `LF{T}`: A LatticeField with the Gaussian distribution.
    """
function lfGaussian(T::DataType, sz::NTuple{N,Integer}, sigma::Number; L::Union{Nothing,Lattice{N}}=nothing, flambda::Real=1) where {N}
    if isnothing(L)
        L = natlat(sz)
    end
    return LF{T}([exp(-(sum(L[i][I[i]]^2 for i = 1:N) / (2 * sigma^2))) for I in CartesianIndices(length.(L))], L, flambda)
end

"""
    lfRing(T::DataType, sz::NTuple{N,Integer}, r::Number, w::Number; L::Union{Nothing, Lattice{N}}=nothing, flambda::Real=1) where N

    Create a LatticeField representing a ring distribution.

    This function generates a ring-shaped intensity distribution, defined on a lattice.

    # Arguments
    - `T::DataType`: The data type of the LatticeField.
    - `sz::NTuple{N,Integer}`: The size of the ring distribution.
    - `r::Number`: The radius of the ring.
    - `w::Number`: The width of the ring.
    - `L::Union{Nothing, Lattice{N}}`: The lattice structure (default: natural lattice).
    - `flambda::Real`: Scaling factor for the wavelength (default: 1).

    # Returns
    - `LF{T}`: A LatticeField with the ring distribution.
    """
function lfRing(T::DataType, sz::NTuple{N,Integer}, r::Number, w::Number; L::Union{Nothing,Lattice{N}}=nothing, flambda::Real=1) where {N}
    if isnothing(L)
        L = natlat(sz)
    end
    return LF{T}([exp(-(sqrt(sum(L[i][I[i]]^2 for i = 1:N)) - r)^2 / (2 * w^2)) for I in CartesianIndices(length.(L))], L, flambda)
end



end