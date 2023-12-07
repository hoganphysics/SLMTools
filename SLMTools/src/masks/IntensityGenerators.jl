module IntensityGenerators
    
using ..LatticeTools: natlat
using FreeTypeAbstraction: findfont, renderface
using Images: Gray, imfilter, Kernel

export makeGaussian, makeRing, makeLetter, makeLetterPattern

"""
    makeGaussian(sz::NTuple{N,Integer}, sigma::Number) where {N}

    Create a multi-dimensional Gaussian distribution.

    # Arguments
    - `sz::NTuple{N,Integer}`: A tuple representing the size of each dimension. Used to determine the 'natural' lattice for the domain. See `natlat` for more information.
    - `sigma::Number`: The standard deviation of the Gaussian. In units of the natural lattice corresponding to `sz`.

    # Returns
    - An N-dimensional array representing a Gaussian distribution.

    This function generates an N-dimensional Gaussian distribution. The distribution is centered in the middle of the array and uses a 'natural' lattice for the domain, scaling the coordinates by 1/sqrt(n) for normalization, which is common in FFT applications.
    """
function makeGaussian(sz::NTuple{N,Integer}, sigma::Number) where {N}
    L = natlat(sz)
    return [exp(-(sum(L[i][I[i]]^2 for i = 1:N) / (2 * sigma^2))) for I in CartesianIndices(length.(L))]
end

"""
    makeRing(sz::NTuple{N,Integer}, r::Number, w::Number) where {N}

    Create a multi-dimensional ring-shaped distribution.

    # Arguments
    - `sz::NTuple{N,Integer}`: A tuple representing the size of each dimension. Used to determine the 'natural' lattice for the domain. See `natlat` for more information.
    - `r::Number`: The radius of the ring. in units of the natural lattice corresponding to `sz`.
    - `w::Number`: The width (standard deviation) of the ring's Gaussian profile. in units of the natural lattice corresponding to `sz`.

    # Returns
    - An N-dimensional array representing a ring-shaped distribution.

    This function creates an N-dimensional distribution shaped like a ring. The ring is defined by a Gaussian profile with a specified radius and width. The distribution is centered in the array and uses a 'natural' lattice for spatial coordinates, scaling by 1/sqrt(n), typical in FFT contexts.
    """
function makeRing(sz::NTuple{N,Integer}, r::Number, w::Number) where {N}
    L = natlat(sz)
    return [exp(-(sqrt(sum(L[i][I[i]]^2 for i = 1:N)) - r)^2 / (2 * w^2)) for I in CartesianIndices(length.(L))]
end

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




end # module IntensityGenerators