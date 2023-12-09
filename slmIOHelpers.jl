# slmIOHelpers.jl
# File input and output functions for working with SLMs and cameras. 
# Requires: Images, FileIO, LatticeFields.jl

#------------------ Image importing and conversion ------------------------------

function getImagesAndFilenames(dir::String,extension::String)
    L = length(extension)                                   # Extension length
    f = s-> length(s) > L && s[end-L+1:end] == extension    # Filter function
    fileNames = filter(f,readdir(dir))
    images = load.(dir .* fileNames)
    return images,fileNames
end

imageToFloatArray(img::Matrix) = convert.(Float64,Gray.(img))             # Converts an image to a matrix of floats, with values in [0,1]. 
itfa = imageToFloatArray;

function castImage(T::DataType,img::Matrix{<:Colorant},L::Lattice{2},flambda::Real)
    return LF{T}(itfa(img),L,flambda)
end
function loadDir(dir::String,extension::String; 
        T::DataType=Intensity, outType::DataType=Float64, L::Union{Lattice{2},Nothing,Real,Tuple{Real,Real}}=nothing, flambda::Real=1.0, 
        cue::Union{String,Char,Nothing}=nothing,look::Union{Symbol,String,Char}=:after)
    # Loads a directory of images, making them into LF{outType} structures, and parses their filenames.  Returns [images], [parameters from filenames].
    images,fileNames = getImagesAndFilenames(dir,extension)
    sz = size(images[1])
    all(size(i) == sz for i in images) || error("Inconsistent image size.")
    if L isa Union{Real,Tuple{Real,Real}}   # Interpret L as the pixel spacing, possibly different in different directions in the case it's a Tuple.
        L = Tuple(1:s for s in sz) .* L
    elseif isnothing(L)
        L = Tuple(1:s for s in sz)
    end
    fields = [castImage(T,i,L,flambda) for i in images]
    if isnothing(cue)
        params = [parseFileName(f; outType=outType) for f in fileNames]
    else
        params = [parseFileName(f,cue,look;outType=outType) for f in fileNames]
    end
    return fields, params
end

#----------------- File name parsing ---------------------------------------------

function parseStringToNum(s::String;outType=nothing)
    # Converts a numeric string to a number.  Numeric strings can only contain numerals and commas. 
        # By default, if there is a comma, then the string is converted to a float.  Otherwise it 
        # is converted to an integer. 
    if isnothing(outType)
        if ',' in s
            s = replace(s,','=>'.')
            return parse(Float64,s)
        else
            return parse(Int,s)
        end
    else
        s = replace(s,','=>'.')
        return parse(outType,s)
    end
end

function parseFileName(name::String,cue::Union{String,Char},look::Union{Symbol,String,Char}=:after; outType=nothing)
    # Finds a numeric string in a filename and turns it into a number.  The numeric string must either precede or follow a 
        # substring or character called "cue".  "look" indicates whether the string follows or precedes the cue.
        # look may be any of :before,:b,"before","b",'b',:after,:a,"after","a",'a'.  All the values starting with a indicate 
        # that the string is after the cue. 
        # The numeric string can only contain numerals and possibly a comma, which stands in for a decimal point.
        # By default, if there's a comma the numeric string is converted to a float, and otherwise it is converted
        # to an integer.  
    look in [:before,:b,"before","b",'b',:after,:a,"after","a",'a'] || throw(ValueError("Unrecognized look value."))
    look in [:before,:b,"before","b",'b'] ? (look = :b) : (look = :a)
    ext = findlast('.',name)
    name = '.'*name[1:ext-1]*'.'    # The periods are to simplify finding the beginning/end of the numeric string. 
    idx = findlast(cue,name)
    if look == :a
        j = findfirst( (x->!(x in [" 0123456789"...,','])) ,name[idx[end]+1:end]) + idx[end]
        s = name[idx[end]+1:j-1]
        return parseStringToNum(s,outType=outType)
    else
        j = findlast( (x->!(x in [" 0123456789"...,','])) ,name[1:idx[1]-1])
        s = name[j+1:idx[1]-1]
        return parseStringToNum(s,outType=outType)
    end
end

function parseFileName(name::String; outType=nothing)
    # Parses a filename to a number.  Assumes that the filename consists of a numeric string followed by a file extension.
    idx = findfirst('.',name)
    return parseStringToNum(name[1:idx-1],outType=outType)
end

#------------------------- LF camera orientation -------------------------------------

function linearFit(xs::Vector{<:Number},ys::Vector{<:Number})
    # Fits a line to the points (x,y) with x in xs, y in ys, via linear regression.  Returns (slope, intercept).
    ((hcat(xs, ones(length(xs))) \ ys)...,)
end
function collapse(x::AbstractArray{T,N},i::Int) where {T<:Number,N}
    # Sums all dimensions of x except the i-th.  
    return sum(x,dims=(1:i-1...,i+1:N...))[:]
end
function clip(x::Number,threshold::Number)
    # Clips a number to zero if it's below the threshold. 
    x>threshold ? x : zero(x)
end
function centroid(img::LF{Intensity,T,N},threshold::Real=0.1) where {T,N}
    # Centroid of a LatticeField, in its own coordinates. 
    clipped = clip.(img.data,threshold)
    s = sum(clipped)
    (s > 0) || error("Black image.  Can't normalize")
    return [(sum(collapse(clipped,i) .* img.L[i]) for i=1:N)...] ./ s
end
function getOrientation(linImgs::Vector{LF{Intensity,T,2}},idxs::Vector{<:Number};roi=nothing, threshold::Number=0.1) where T<:Number
    if !isnothing(roi)
        linImgs = [i[roi...] for i in linImgs]
    end
    cs = vcat( (centroid(i,threshold)' for i in linImgs )... );
    xs, xi, ys, yi = linearFit(idxs,cs[:,1])...,linearFit(idxs,cs[:,2])...
    theta = angle(xs + im*ys);
    return [xi,yi],theta
end

#--------------------------- Dualization --------------------------------------------

function dualate(f::LF{S,T,2},L::Lattice{2},center::Vector{<:Real},θ::Real,flambda::Real=1.0; 
        roi=nothing, interpolation=cubic_spline_interpolation, naturalize=false,bc=zero(T)) where {S<:FieldVal,T<:Number}
    # Corrects for offset and tilt of f.L and then interpolates f onto a lattice dual to L. 
    length(center) == length(f.L) || error("Incompatible lengths for center and f.L.")
    if !isnothing(roi)
        f = f[roi...]
    end
    L0 = Tuple(f.L[i] .- center[i] for i=1:2)
    interp = interpolation(L0, f.data, extrapolation_bc=bc)
    dL = dualShiftLattice(L,flambda)
    Lx0, Lxs, Ly0, Lys = dL[1][1], step(dL[1]), dL[2][1], step(dL[2])
    R = [cos(θ) -sin(θ) ; sin(θ) cos(θ)]
    r0, dx, dy = (R * [Lx0, Ly0]), Lxs * R[:,1], Lys * R[:,2]
    if naturalize
        return LF{S}([interp((r0 + dx*I[1] + dy*I[2])...) for I in CartesianIndices(length.(dL))],natlat(length.(dL)),1.0);
    else
        return LF{S}([interp((r0 + dx*I[1] + dy*I[2])...) for I in CartesianIndices(length.(dL))],dL,flambda);
    end
end
dualate(fs::Vector{<:LF{S}},L::Lattice{2},center::Vector{<:Real},theta::Real,flambda::Real=1.0;kwargs...) where {S<:FieldVal} = 
    [dualate(f,L,center,theta,flambda;kwargs...) for f in fs]

#--------------------------- File output --------------------------------------------

function savePhase(x::AbstractArray{ComplexF64,2},name::String)
    save(name,Gray.((angle.(x) .+ pi) ./ (2pi)))
end
function savePhase(x::LF{ComplexPhase},name::String)
    savePhase(x.data,name)
end
function savePhase(x::LF{RealPhase},name::String)
    save(name,Gray.(mod.(x.data,1)))
end

#------------------------- Viewing in a notebook ------------------------------------
cycle1(x::Complex) = (angle(x) + pi)/(2pi)
function look(f::LF{RealPhase})
    return Gray.((f.data .- minimum(f.data)) .% 1)
end
function look(f::LF{ComplexPhase})
    return Gray.(cycle1.(f.data))
end
function look(f::LF{Modulus})
    return Gray.(f.data ./ maximum(f.data))
end
function look(f::LF{ComplexAmp})
    return hcat(Gray.(abs.(f.data) ./ maximum(abs.(f.data))),Gray.(cycle1.(f.data)))
end
function look(f::LF{Intensity})
    return Gray.(f.data ./ maximum(f.data))
end

function look(x::LF...)
    return hcat((look(y) for y in x)...)
end

#------------------- Shifted Fourier transforms and natlats -------------------------------

sft(v) = fftshift(fft(ifftshift(v)))
isft(v) = fftshift(ifft(ifftshift(v)))
sft(v::LF{ComplexAmplitude}) = LF{ComplexAmplitude}(fftshift(fft(ifftshift(v.data))),dualShiftLattice(v.L,v.flambda),v.flambda)
isft(v::LF{ComplexAmplitude}) = LF{ComplexAmplitude}(fftshift(ifft(ifftshift(v.data))),dualShiftLattice(v.L,v.flambda),v.flambda)

function natlat(n::Int)
    return (-floor(n/2):floor((n-1)/2)) ./ sqrt(n)
end
function natlat(s::NTuple{N,Integer}) where N
    return ((natlat(n) for n in s)...,)
end;

function naturalize(f::LF{S,T,N}) where {S<:FieldVal,T,N}
    # Converts LatticeField to be difined on the natural lattice of the same size. 
    return LF{S}(f.data,natlat(size(f)),1.0)
end

#------------------ Standard phases and intensities --------------------------------------


function lfRampedParabola(T::DataType,imSize::NTuple{N,Integer}, pixelScale::Real, center::NTuple{N,Real}, pAmp::Real, rAmps::NTuple{N,Real};L::Union{Nothing,Lattice{N}}=nothing,flambda::Real=1) where N
    if isnothing(L)
        L = natlat(imSize)
    end
    im = ([ pAmp * sum( (Tuple(I) .- center).^2 ) / (2*pixelScale^2) for I in CartesianIndices(imSize)] .+ 
          [ sum( (Tuple(I) .- center) .* rAmps ) / pixelScale  for I in CartesianIndices(imSize)]) |> (x->mod.(x,1))
    return LF{T}(im,L,flambda)
end


function lfGaussian(T::DataType,sz::NTuple{N,Integer},sigma::Number;L::Union{Nothing,Lattice{N}}=nothing,flambda::Real=1) where N
    if isnothing(L)
        L = natlat(sz)
    end
    return LF{T}([exp(-(sum(L[i][I[i]]^2 for i=1:N)/(2*sigma^2))) for I in CartesianIndices(length.(L))],L,flambda)
end

function lfRing(T::DataType,sz::NTuple{N,Integer},r::Number,w::Number;L::Union{Nothing,Lattice{N}}=nothing,flambda::Real=1) where N
    if isnothing(L)
        L = natlat(sz)
    end
    return LF{T}([exp(-(sqrt(sum(L[i][I[i]]^2 for i=1:N))-r)^2/(2*w^2)) for I in CartesianIndices(length.(L))],L,flambda)
end

