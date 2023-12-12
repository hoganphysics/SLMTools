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

# TODO these are older functions, need to be updated to use LatticeFields
"""
    makeLetter(c::Char, sz::Tuple{Int,Int}, fnt="arial bold"; returnImg=true)

    Create an image of a specified character in a given font.

    # Arguments
    - `c::Char`: The character to render.
    - `sz::Tuple{Int,Int}`: The size of the image (rows, columns).
    - `fnt`: The font used to render the character, defaults to "arial bold".

    # Keyword Arguments
    - `returnImg`: If `true`, returns the image as a `Gray` array; if `false`, returns the raw numerical array.

    # Returns
    - An image array of the specified character.

    This function renders a single character `c` using the specified font `fnt` and creates an image of size `sz`. The rendered character is centered in the image. The image is returned as a grayscale image if `returnImg` is true, otherwise as a raw numerical array.
    """
function makeLetter(c::Char, sz::Tuple{Int,Int}, fnt="arial bold"; returnImg=true)
    f = findfont(fnt)
    l = renderface(f, c, sz[1])[1]
    img = zeros(Float64, sz)
    img[CartesianIndices(size(l)).+CartesianIndex((sz .- size(l)) .÷ 2)] .= convert.(Float64, l) ./ 255
    returnImg ? (return Gray.(img)') : (return img')
end

"""
    makeLetterPattern(tw::Int, NT::Int, word::String)

    Create an image of a word using a tiled pattern of letters.

    # Arguments
    - `tw::Int`: The size of each letter in the pattern.
    - `NT::Int`: The size of the image (rows, columns).
    - `word::String`: The word to render.

    # Returns
    - An image array of the specified word.

    This function renders a word using a tiled pattern of letters. The letters are rendered using `makeLetter` and are tiled in a grid pattern to form the word. The image is returned as a grayscale image.
    """
function makeLetterPattern(tw::Int, NT::Int, word::String)
    txt = hcat([makeLetter(i, (tw, tw))[:, (i == 'I' ? (floor(Int, 0.35 * tw):floor(Int, 0.65 * tw)) : (floor(Int, 0.2 * tw):floor(Int, 0.85 * tw)))] for i in word]...)
    arr = zeros(typeof(txt[1]), NT, NT)
    arr[CartesianIndices(size(txt)).+CartesianIndex((size(arr) .- size(txt)) .÷ 2).+0*CartesianIndex(NT ÷ 4, 0)] .= txt
    arr = imfilter(arr, Kernel.gaussian(10))
    return arr
end



end