module LFTemplates

using ..LatticeTools
using FreeTypeAbstraction: findfont, renderstring!

export lfRampedParabola, lfGaussian, lfRing, lfParabolaCap, ftaText, lfText, lfRect, lfRand


################################## Helper functions #################################

function lfStandardOutputFormat(T::DataType,data::Array{S,N},L::Lattice{N},flambda::Real) where {S,N}
    if T==ComplexPhase
        if S<:Real
            return wrap(LF{RealPhase}(data,L,flambda))
        elseif S<:Complex
            return LF{ComplexPhase}(data,L,flambda)
        else
            error("`data` type not understood.")
        end
    elseif T in (Intensity,Modulus)
        if any(data .< 0)
            println("Warning: negative values in nominally non-negative LF data field. Clipping to zero.")
        end
        return LF{T}(ramp.(data),L,flambda)
    else
        return LF{T}(data,L,flambda)
    end
end

l2form(L::Lattice{N},M::Matrix{T}) where {N,T<:Number} = .+( (M[i,j] .* toDim(L[i],i,N) .* toDim(L[j],j,N) for i=1:N,j=1:N)... )
ramp(x::T) where {T<:Number} = (x < 0 ? zero(T) : x)

"""
    ftaText(str::String,sz::Tuple{Int,Int}; fnt = "arial bold",pixelsize::Union{Int,Nothing}=nothing,halign=:hcenter,valign=:vcenter,options...)

    Make text string `str` into a float array of size `sz`.  The letter size is set with optional kwarg pixelsize. 
    WARNING: This function (and the package it derives from) does not seem to work on Linux. 

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


################################## Metaprogramming for auto-generating methods #################################

methodPattern1 = quote
    function fname(T::DataType,L::Lattice{N} ; 
            center::NTuple{N,Real}=Tuple(0.0 for i=1:N), flambda::Real=1.0) where {N}
        Lc = Tuple(L[i] .- center[i] for i=1:N)
        begin
        end
    end
end

methodPattern2 = quote
    function fname(lfPattern::LF{T,S,N} ; center::NTuple{N,Real}=Tuple(0.0 for i=1:N)) where {T<:FieldVal,S,N}
        L = lfPattern.L
        flambda = lfPattern.flambda
        Lc = Tuple(L[i] .- center[i] for i=1:N)
        begin
        end
    end
end

patternRef_args(x::Expr)         = x.args[2].args[1].args[1].args     # push! arguments here
patternRef_where(x::Expr)        = x.args[2].args[1].args             # push! type parameters here
patternRef_body(x::Expr)         = x.args[2].args[2].args             # Update [end] of this to function body

extract_func(x::Expr) = (x.args[1].head == :call) ? (x.args[1]) : (x.args[1].args[1])
funcRef_name(x::Expr) = extract_func(x).args[1]
funcRef_posArgs(x::Expr) = extract_func(x).args[2].head == :parameters ? extract_func(x).args[3:end] : extract_func(x).args[2:end]
funcRef_kwArgs(x::Expr) = extract_func(x).args[2].head == :parameters ? extract_func(x).args[2].args : []
funcRef_params(x::Expr) = (x.args[1].head == :call) ? [] : (x.args[1].args[2:end])
funcRef_body(x::Expr) = x.args[2]

function updateFuncExpr!(pattern,func)
    patternRef_args(pattern)[1] = funcRef_name(func)                  # Set function name
    push!(patternRef_args(pattern),funcRef_posArgs(func)...)          # Set positional args
    push!(patternRef_args(pattern)[2].args,funcRef_kwArgs(func)...)   # Set kwargs
    push!(patternRef_where(pattern),funcRef_params(func)...)          # Set parametric type info
    patternRef_body(pattern)[end] = funcRef_body(func)                # Set function body
    push!(patternRef_body(pattern),:(return lfStandardOutputFormat(T,p,L,flambda)))  # Generic return statement
end

macro addTemplateMethods(funcExpr)
    patt1 = copy(methodPattern1)
    updateFuncExpr!(patt1,funcExpr)
    
    patt2 = copy(methodPattern2)
    updateFuncExpr!(patt2,funcExpr)

    return quote 
        $(esc(patt1))
        $(esc(patt2))
    end
end

################################# LF Templates #################################

@addTemplateMethods function lfRand(; R::DataType=Float64)
    p = rand(R,length.(L))
end

@addTemplateMethods function lfParabola(quad::Real, lin::NTuple{N,Real}=Tuple(0.0 for i=1:N))
    p = quad/2 * r2(Lc) + ldot(lin,Lc)
end

@addTemplateMethods function lfParabola(quad::Matrix{R}, lin::NTuple{N,Real}=Tuple(0.0 for i=1:N)) where {R<:Real}
    p = l2form(Lc,quad) ./ 2 + ldot(lin,Lc)
end

@addTemplateMethods function lfGaussian(radius::Real, norm::Real=1.0)
    p = exp.(-r2(Lc) ./(2*radius^2))
    p .*= sqrt(norm / (sum(p.^2) * *(step.(Lc)...)))
end

@addTemplateMethods function lfGaussian(covar::Matrix{R}, norm::Real=1.0) where {R}
    p = exp.(-l2form(Lc,covar) ./ 2)
    p .*= sqrt(norm / (sum(p.^2) * *(step.(Lc)...)))
end

@addTemplateMethods function lfRing(radius::Number, width::Number)
    p = exp.( -(sqrt.(r2(Lc)) .- radius).^2 ./ (2*width^2) )
end

@addTemplateMethods function lfCap(curvature::Real, height::Real)
    p = ramp.(height .- curvature*r2(L)/2)
end

@addTemplateMethods function lfRect(sides::NTuple{N,Real},height::Real=1.0)
	p = zeros(Float64,length.(Lc))
	p[(abs.(Lc[i]) .<= sides[i]/2 + eps() for i=1:N)...] .= height
end

@addTemplateMethods function lfText(str::String; R::DataType=Float64, pixelsize::Union{Int,Nothing}=nothing, 
        fnt = "arial bold", halign=:hcenter, valign=:vcenter, options...)
    # WARNING: This function probably doesn't work on Linux machines, due to a bug in the FreeTypeAbstraction package.
    p = convert.(R,ftaText(str,length.(L);pixelsize=pixelsize,fnt=fnt,halign=halign,valign=valign, options...))
end


end
