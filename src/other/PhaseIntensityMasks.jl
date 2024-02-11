module PhaseIntensityMasks

using ..LatticeTools
using FreeTypeAbstraction: findfont, renderstring!

export lfRampedParabola, lfGaussian, lfRing, lfParabolaCap, ftaText, lfText, lfRect, lfRand

"""
    lfRampedParabola(T::DataType, imSize::NTuple{N,Integer}, pixelScale::Real, center::NTuple{N,Real}, 
                     pAmp::Real, rAmps::NTuple{N,Real}; L::Union{Nothing, Lattice{N}}=nothing, flambda::Real=1) where N

    Create a LatticeField representing a ramped parabolic distribution.  The distribution is 
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
    lfRampedParabola(T::DataType, L::Lattice{N}, pAmp::Real, rAmps::NTuple{N,Real}, offset::Real=0.0; lambda::Real=1.0) where N

    Create a LatticeField representing a ramped parabolic distribution.  Accepts input lattice and polynomial coefficients. 
	
    # Arguments
    - `T::DataType`: The data type of the LatticeField.
    - `L::Lattice{N}`: Lattice for the LF. 
    - `pAmp::Real`: The amplitude of the parabola.
    - `rAmps::NTuple{N,Real}`: The linear coefficients in each dimension.
    - `flambda::Real`: Scaling factor for the wavelength (default: 1).

    # Returns
    - `LF{T}`: A LatticeField with the ramped parabolic distribution.
    """
function lfRampedParabola(T::DataType, L::Lattice{N}, pAmp::Real, rAmps::NTuple{N,Real}, offset::Real=0.0; flambda::Real=1.0) where N
    im = pAmp*r2(L)/2 .+ ldot(rAmps,L) .+ offset
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

ramp(x::T) where {T<:Number} = (x < 0 ? zero(T) : x)

"""
    lfParabolaCap(T::DataType, L::Lattice{N}, curvature::Real, height::Real, flambda::Real=1.0; center::Union{Nothing,Vector{<:Real}}=nothing) where N

    Create a LatticeField representing a parabolic cap, i.e. an inverted parabola truncated at zero.  
	The formula is `ramp.(height .- curvature*r2(L)/2)`.  
	Can be recentered to any location via the optional `center` kwarg. 

    # Arguments
    - `T::DataType`: The data type of the LatticeField.
    - `L::Lattice{N}`: Lattice.
    - `curvature::Real`: The quadratic coefficient of the parabola, defined such that the parabola is `curvature*x^2/2`. 
    - `height::Real`: Cap height. 
    - `flambda::Real=1.0`: Scale factor.
    - `center::Union{Nothing,Vector{<:Real}}=nothing`: Center of the cap.

    # Returns
    - `LF{T}`: A LatticeField with the parabolic cap distribution.
    """
function lfParabolaCap(T::DataType, L::Lattice{N}, curvature::Real, height::Real, flambda::Real=1.0; center::Union{Nothing,Vector{<:Real}}=nothing) where N
    if isnothing(center)
        center = zeros(N)
    end
    L = Tuple(L[i] .- center[i] for i=1:N)
    return LF{T}(ramp.(height .- curvature*r2(L)/2),L,flambda)
end

"""
    ftaText(str::String,sz::Tuple{Int,Int}; fnt = "arial bold",pixelsize::Union{Int,Nothing}=nothing,halign=:hcenter,valign=:vcenter,options...)

    Make text string `str` into a float array of size `sz`.  The letter size is set with optional kwarg pixelsize. 

    # Arguments
    - `str::String`: Text string.
    - `sz::NTuple{N,Integer}`: Desired array size.
    - `fnt = "arial bold"`: Font.
    - `pixelsize::Union{Int,Nothing}=nothing`: Text size parameter.
    - Various other options: You should know what you're doing to set these. 

    # Returns
    - An array displaying the text.
    """
function ftaText(str::String,sz::Tuple{Int,Int}; fnt = "arial bold",pixelsize::Union{Int,Nothing}=nothing,halign=:hcenter,valign=:vcenter,options...)
    if isnothing(pixelsize)
        pixelsize = sz[2] ÷ length(str)
    end
    face = findfont(fnt)
    x0, y0 = sz .÷ 2
    arr = zeros(UInt8,sz...)    # Text will go here
    renderstring!(arr,str,face,pixelsize, x0, y0; halign=halign, valign=valign, options...)
    return convert.(Float64,arr)./255
end

"""
    lfText(S::DataType,str::String,sz::Union{Nothing,Tuple{Int,Int}}=nothing; T::DataType=Float64,
        L::Union{Lattice,Nothing}=nothing, flambda::Union{Real,Nothing}=nothing, lfTemplate::Union{LF,Nothing}=nothing,
        pixelsize::Union{Int,Nothing}=nothing, fnt = "arial bold", halign=:hcenter, valign=:vcenter, options...)

    Make text string `str` into LatticeField of type `S`.  The letter size is set with optional kwarg pixelsize. 
	You must set the output size either via the positional arg sz, the kwarg L (setting the lattice explicitly),
	or the kwarg lfTemplate, which provides an existing LF from which to get the lattice and flambda. 

    # Arguments
    - `str::String`: Text string.
    - `sz::NTuple{N,Integer}`: Desired array size.
    - `fnt = "arial bold"`: Font.
    - `pixelsize::Union{Int,Nothing}=nothing`: Text size parameter.
	- `L`: A lattice.
	- `flambda`: Scale factor.
	- `lfTemplate`: A LF from which to get the lattice and flambda.  This overrides L and flambda kwargs. 
    - Various other options: You should know what you're doing to set these. 

    # Returns
    - An array displaying the text.
    """
function lfText(S::DataType,str::String,sz::Union{Nothing,Tuple{Int,Int}}=nothing; T::DataType=Float64,
        L::Union{Lattice,Nothing}=nothing, flambda::Union{Real,Nothing}=nothing, lfTemplate::Union{LF,Nothing}=nothing,
        pixelsize::Union{Int,Nothing}=nothing, fnt = "arial bold", halign=:hcenter, valign=:vcenter, options...)
    if !isnothing(lfTemplate)
        flambda = lfTemplate.flambda
        L = lfTemplate.L
    end
    if isnothing(sz) && isnothing(L)
        error("Must specify text size either via sz or a lattice L")
    end
    if isnothing(sz)
        sz = length.(L)
    end
    if isnothing(L)
        L = natlat(sz)
    end
    if isnothing(flambda)
        flambda = 1
    end
    txt = convert.(T,ftaText(str,sz;pixelsize=pixelsize,fnt=fnt,halign=halign,valign=valign, options...))
    return LF{S}(txt,L,flambda)
end

"""
    lfRect(lft::LF{S,T,N},s::NTuple{N,Real},c::NTuple{N,Real}=Tuple(zeros(N))) where {S,T,N}

    Make LatticeField in the shape of a rectangle.  `s` specifies the side length in each direction, and
	`c` specifies the center.  `lft` is a template from which to infer the lattice and flambda. 

    # Arguments
    - `lft`: Template LatticeField.
    - `s::NTuple{N,Real}`: Rectangle side lengths.
    - `c::NTuple{N,Real}`: Rectangle center.

    # Returns
    - A LatticeField{S} with rectangular field profile.
    """
function lfRect(lft::LF{S,T,N},s::NTuple{N,Real},c::NTuple{N,Real}=Tuple(zeros(N))) where {S,T,N}
	x=zeros(T,size(lft))
	x[(abs.(lft.L[i].-c[i]) .< s[i] for i=1:N)...] .= 1
    return LF{S}( x, lft.L,lft.flambda )
end


"""
    lfRand(S::DataType,T::DataType,L::Lattice{N},flambda::Real) where N

    Make LatticeField with random data.

    # Arguments
    - `S::DataType`: LatticeField FieldVal.
    - `T::DataType`: Data type (i.e. Float64, ComplexF64, etc.).
    - `L::Lattice{N}`: Lattice.
	- `flambda::Real`: flambda value. 

    # Returns
    - A LatticeField{S} with random field profile.
	"""
lfRand(S::DataType,T::DataType,L::Lattice{N},flambda::Real) where N = LF{S}(rand(T,length.(L)),L,flambda)

"""
    lfRand(lft::LatticeField{S,T,N}) where {S,T,N}

    Make LatticeField with random data from LF template.  The output inherits FieldVal, lattice, and flambda
	from the template. 

    # Arguments
    - `lft::LatticeField{S,T,N}`: LatticeField template.

    # Returns
    - A LatticeField{S,T,N} with random field profile.
	"""
lfRand(lft::LatticeField{S,T,N}) where {S,T,N} = LF{S}(rand(T,size(lft)),lft.L,lft.flambda)

#--------------- Deprecated text functions --------------------------------------------

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