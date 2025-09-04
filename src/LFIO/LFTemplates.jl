#=
Template functions for generating common LatticeField distributions.  These are implemented with metaprogramming to reduce boilerplate.

Requires:
    using FreeTypeAbstraction: findfont, renderstring!
=#

export lfParabola, lfGaussian, lfRing, lfCap, ftaText, lfText, lfRect, lfRand, lfHeart, lfSmile, lfPointer, lfBlur

################################## Function docstrings #################################

"""
    lfParabola(...,quad::Real, lin::NTuple{N,Real}=Tuple(0.0 for i=1:N))

LF template function generating a parabolic LatticeField with specified quadratic and linear coefficients. 

All LF template functions have a method which accepts a Lattice and LF type (e.g. Intensity, Modulus) and another
method which acccepts an LF from which the Lattice and LF type are inferred. In either case, the first argument
or arguments are these template objects, and subsequent arguments are peculiar to the function.  See examples 
below for how to use the template arguments.  In this case, the peculiar arguments are as follows: 

# Peculiar arguments
- `quad::Real`: Second derivative of the parabola.
- `lin::NTuple{N,Real}`: First derivative of the parabola in each dimension.

# Template arguments
Use one or the other of the following argument combinations as the initial argument(s) to this function:
- `T::FieldVal`, `L::Lattice`: The fieldval (e.g. `RealPhase``, `Intensity`) and lattice.
- `f::LF{T}`: An LF whose type and lattice fill the roles of `T`, `L` in the above version. 

# Keyword arguments
- `center::NTuple{N,Real}=Tuple(0.0 for i=1:N)`: Center position by which to offset the output. 
- `flambda::Real=1.0`: flambda value.  This is only available when the template arguments are of the `T,L` form. 

# Returns
- `LatticeField` with a parabolic distribution.

# Notes
- The output is calculated as `quad * r2(L)/2 + ldot(lin,L)`, where `r2(L)` is the squared radius
    in the lattice coordinates and `ldot` is the lattice dot product.

# Examples
- f1 = lfParabola(Intensity,natlat(128,128),1.5)
- f2 = lfParabola(f1,2.1,(0.1,0.2))
"""
function lfParabola end

"""
    lfGaussian(...,radius::Real, norm::Real=1.0)

LF template function generating a Gaussian LatticeField with specified radius and normalization. 

All LF template functions have a method which accepts a Lattice and LF type (e.g. Intensity, Modulus) and another
method which acccepts an LF from which the Lattice and LF type are inferred. In either case, the first argument
or arguments are these template objects, and subsequent arguments are peculiar to the function.  See examples 
below for how to use the template arguments.  In this case, the peculiar arguments are as follows: 

# Peculiar arguments
- `radius::Real`: Standard deviation of the Gaussian.
- `norm::Real`: Normalization factor for the Gaussian.

# Template arguments
Use one or the other of the following argument combinations as the initial argument(s) to this function:
- `T::FieldVal`, `L::Lattice`: The fieldval (e.g. `RealPhase``, `Intensity`) and lattice.
- `f::LF{T}`: An LF whose type and lattice fill the roles of `T`, `L` in the above version. 

# Keyword arguments
- `center::NTuple{N,Real}=Tuple(0.0 for i=1:N)`: Center position by which to offset the output. 
- `flambda::Real=1.0`: flambda value.  This is only available when the template arguments are of the `T,L` form. 

# Returns
- `LatticeField` with a Gaussian distribution.

# Notes
- The output is calculated as `exp(-r2(L)/(2*radius^2))`, where `r2(L)` is the squared radius
    in the lattice coordinates.

# Examples
- f1 = lfGaussian(Intensity,natlat(128,128),1.5)
- f2 = lfGaussian(f1,2.1)
"""
function lfGaussian end

"""
    lfRing(<template args>,radius::Number, width::Number; <kwargs>)

LF template function generating a ring-shaped LatticeField with specified radius and width.  
    
All LF template functions have a method which accepts a Lattice and LF type (e.g. Intensity, Modulus) and another
method which acccepts an LF from which the Lattice and LF type are inferred. In either case, the first argument
or arguments are these template objects, and subsequent arguments are peculiar to the function.  See examples 
below for how to use the template arguments.  In this case, the peculiar arguments are as follows: 

# Peculiar arguments
- `radius::Number`: Radius of the ring.
- `width::Number`: Width of the ring.

# Template arguments
Use one or the other of the following argument combinations as the initial argument(s) to this function:
- `T::FieldVal`, `L::Lattice`: The fieldval (e.g. `RealPhase``, `Intensity`) and lattice.
- `f::LF{T}`: An LF whose type and lattice fill the roles of `T`, `L` in the above version. 

# Keyword arguments
- `center::NTuple{N,Real}=Tuple(0.0 for i=1:N)`: Center position by which to offset the output. 
- `flambda::Real=1.0`: flambda value.  This is only available when the template arguments are of the `T,L` form. 

# Returns
- `LatticeField` with a ring-shaped distribution.

# Notes
- The output is calculated as `exp(-(sqrt(r2(L)) - radius)^2 / (2*width^2))`, where `r2(L)` is the squared radius
    in the lattice coordinates.

# Examples
- f1 = lfRing(Intensity,natlat(128,128),2.0,0.5)
- f2 = lfRing(f1,2.0,0.5)
"""
function lfRing end

"""
    lfCap(...,curvature::Real, height::Real)

LF template function generating a truncated parabolic LatticeField with specified curvature and height.  

All LF template functions have a method which accepts a Lattice and LF type (e.g. Intensity, Modulus) and another
method which acccepts an LF from which the Lattice and LF type are inferred. In either case, the first argument
or arguments are these template objects, and subsequent arguments are peculiar to the function.  See examples 
below for how to use the template arguments.  In this case, the peculiar arguments are as follows: 

# Peculiar arguments
- `curvature::Real`: Second derivative of the parabola.
- `height::Real`: Height of the cap.

# Template arguments
Use one or the other of the following argument combinations as the initial argument(s) to this function:
- `T::FieldVal`, `L::Lattice`: The fieldval (e.g. `RealPhase``, `Intensity`) and lattice.
- `f::LF{T}`: An LF whose type and lattice fill the roles of `T`, `L` in the above version. 

# Keyword arguments
- `center::NTuple{N,Real}=Tuple(0.0 for i=1:N)`: Center position by which to offset the output. 
- `flambda::Real=1.0`: flambda value.  This is only available when the template arguments are of the `T,L` form. 

# Returns
- `LatticeField` with a parabola cap distribution.

# Notes
- The output is calculated as `ramp(curvature * r2(L)/2)`, where `r2(L)` is the squared radius
    in the lattice coordinates.

# Examples
- f1 = lfCap(Intensity,natlat(128,128),1.5,2.0)
- f2 = lfCap(f1,2.1,2.0)
"""
function lfCap end

"""
    lfText(...,str::String)

LF template function generating a text LatticeField with specified string. 

All LF template functions have a method which accepts a Lattice and LF type (e.g. Intensity, Modulus) and another
method which acccepts an LF from which the Lattice and LF type are inferred. In either case, the first argument
or arguments are these template objects, and subsequent arguments are peculiar to the function.  See examples 
below for how to use the template arguments.  In this case, the peculiar arguments are as follows: 

# Peculiar arguments
- `str::String`: Text string to be rendered.
- `pixelsize::Union{Int,Nothing}=nothing`: Text size parameter.
- `fnt = "arial bold"`: Font.
- `halign=:hcenter`: Horizontal alignment.
- `valign=:vcenter`: Vertical alignment.
- `options...`: Additional rendering options.

# Template arguments
Use one or the other of the following argument combinations as the initial argument(s) to this function:
- `T::FieldVal`, `L::Lattice`: The fieldval (e.g. `RealPhase``, `Intensity`) and lattice.
- `f::LF{T}`: An LF whose type and lattice fill the roles of `T`, `L` in the above version. 

# Keyword arguments
- `center::NTuple{N,Real}=Tuple(0.0 for i=1:N)`: Center position by which to offset the output. 
- `flambda::Real=1.0`: flambda value.  This is only available when the template arguments are of the `T,L` form. 

# Returns
- `LatticeField` with a text distribution.

# Notes
- The output is generated using the `ftaText` function, which renders the text string into a float array.

# Examples
- f1 = lfText(Intensity,natlat(128,128),"C")
- f2 = lfText(f1,"C")
"""
function lfText end

    """
    lfRect(...,sides::NTuple{N,Real},height::Real=1.0)

LF template function generating a rectangular shaped (i.e. boxcar function) LatticeField.  

All LF template functions have a method which accepts a Lattice and LF type (e.g. Intensity, Modulus) and another
method which acccepts an LF from which the Lattice and LF type are inferred. In either case, the first argument
or arguments are these template objects, and subsequent arguments are peculiar to the function.  See examples 
below for how to use the template arguments.  In this case, the peculiar arguments are as follows: 

# Peculiar arguments
- sides::NTuple{N,Real}: Lengths of the sides of the rectangle in each dimension.
- height::Real=1.0: Height of the boxcar function.

# Template arguments
Use one or the other of the following argument combinations as the initial argument(s) to this function:
- `T::FieldVal`, `L::Lattice`: The fieldval (e.g. `RealPhase``, `Intensity`) and lattice.
- `f::LF{T}`: An LF whose type and lattice fill the roles of `T`, `L` in the above version. 

# Keyword arguments
- `center::NTuple{N,Real}=Tuple(0.0 for i=1:N)`: Center position by which to offset the output. 
- `flambda::Real=1.0`: flambda value.  This is only available when the template arguments are of the `T,L` form. 

# Returns
- `LatticeField` with a rectangular shape.

# Examples
- f1 = lfRect(Intensity,natlat(128,128),(1.5,2.3))
- f2 = lfRect(f1,(1.5,1.3))
"""
function lfRect end

"""
    lfRand(...)

LF template function generating a random LatticeField. 

All LF template functions have a method which accepts a Lattice and LF type (e.g. Intensity, Modulus) and another
method which acccepts an LF from which the Lattice and LF type are inferred. In either case, the first argument
or arguments are these template objects, and subsequent arguments are peculiar to the function.  See examples 
below for how to use the template arguments.  In this case, there are no peculiar arguments, save an optional 
keyword argument to specify the data type of the field values.

# Peculiar arguments
- R::DataType=Float64: Data type of the random values.

# Template arguments
Use one or the other of the following argument combinations as the initial argument(s) to this function:
- `T::FieldVal`, `L::Lattice`: The fieldval (e.g. `RealPhase``, `Intensity`) and lattice.
- `f::LF{T}`: An LF whose type and lattice fill the roles of `T`, `L` in the above version. 

# Keyword arguments
- `center::NTuple{N,Real}=Tuple(0.0 for i=1:N)`: Center position by which to offset the output. 
- `flambda::Real=1.0`: flambda value.  This is only available when the template arguments are of the `T,L` form. 

# Returns
- `LatticeField` with random values.

# Examples
- f1 = lfRand(Intensity,natlat(128,128))
- f2 = lfRand(f1)
"""
function lfRand end

"""
    lfHeart(...,scale; flip=false)

LF template function generating a heart-shaped LatticeField.

All LF template functions have a method which accepts a Lattice and LF type (e.g. Intensity, Modulus) and another
method which acccepts an LF from which the Lattice and LF type are inferred. In either case, the first argument
or arguments are these template objects, and subsequent arguments are peculiar to the function.  See examples 
below for how to use the template arguments.  In this case, the peculiar arguments are as follows:  

# Peculiar arguments
- `scale`: Overall scale of the heart shape.
- `flip=false`: If true, flip the heart vertically.

# Template arguments
Use one or the other of the following argument combinations as the initial argument(s) to this function:
- `T::FieldVal`, `L::Lattice`: The fieldval (e.g. `RealPhase``, `Intensity`) and lattice.
- `f::LF{T}`: An LF whose type and lattice fill the roles of `T`, `L` in the above version.

# Keyword arguments
- `center::NTuple{N,Real}=Tuple(0.0 for i=1:N)`: Center position by which to offset the output. 
- `flambda::Real=1.0`: flambda value.  This is only available when the template arguments are of the `T,L` form.

# Returns
- `LatticeField` with a heart shape.

# Examples
- f1 = lfHeart(Intensity,natlat(128,128),1.5)
- f2 = lfHeart(f1,2.1)
"""
function lfHeart end

"""
    lfSmile(...,scale; flip=false)

LF template function generating a smiley face LatticeField.

All LF template functions have a method which accepts a Lattice and LF type (e.g. Intensity, Modulus) and another
method which acccepts an LF from which the Lattice and LF type are inferred. In either case, the first argument
or arguments are these template objects, and subsequent arguments are peculiar to the function.  See examples
below for how to use the template arguments.  In this case, the peculiar arguments are as follows:

# Peculiar arguments
- `scale`: Overall scale of the smiley face.
- `flip=false`: If true, flip the smiley vertically.

# Template arguments
Use one or the other of the following argument combinations as the initial argument(s) to this function:
- `T::FieldVal`, `L::Lattice`: The fieldval (e.g. `RealPhase``, `Intensity`) and lattice.
- `f::LF{T}`: An LF whose type and lattice fill the roles of `T`, `L` in the above version.

# Keyword arguments
- `center::NTuple{N,Real}=Tuple(0.0 for i=1:N)`: Center position by which to offset the output.
- `flambda::Real=1.0`: flambda value.  This is only available when the template arguments are of the `T,L` form.

# Returns
- `LatticeField` with a smiley face.

# Examples
- f1 = lfSmile(Intensity,natlat(128,128),1.5)
- f2 = lfSmile(f1,2.1)
"""
function lfSmile end

"""
    lfPointer(...,scale; flip=false)

LF template function generating a pointing hand emoji LatticeField.

All LF template functions have a method which accepts a Lattice and LF type (e.g. Intensity, Modulus) and another
method which acccepts an LF from which the Lattice and LF type are inferred. In either case, the first argument
or arguments are these template objects, and subsequent arguments are peculiar to the function.  See examples
below for how to use the template arguments.  In this case, the peculiar arguments are as follows:

# Peculiar arguments
- `scale`: Overall scale of the pointing hand.
- `flip=false`: If true, flip the hand vertically.

# Template arguments
Use one or the other of the following argument combinations as the initial argument(s) to this function:
- `T::FieldVal`, `L::Lattice`: The fieldval (e.g. `RealPhase``, `Intensity`) and lattice.
- `f::LF{T}`: An LF whose type and lattice fill the roles of `T`, `L` in the above version.

# Keyword arguments
- `center::NTuple{N,Real}=Tuple(0.0 for i=1:N)`: Center position by which to offset the output.
- `flambda::Real=1.0`: flambda value.  This is only available when the template arguments are of the `T,L` form.

# Returns
- `LatticeField` with a pointing hand emoji.

# Examples
- f1 = lfPointer(Intensity,natlat(128,128),1.5)
- f2 = lfPointer(f1,2.1)
"""
function lfPointer end

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
# ramp(x::T) where {T<:Number} = (x < 0 ? zero(T) : x)

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

################################## lfBlur #################################

"""
    lfBlur(f::LF{T,S,N},r::Real) where {T,S,N}

Apply a Gaussian blur to the LatticeField `f` with standard deviation `r`.

# Arguments
- `f::LF{T,S,N}`: Input LatticeField.
- `r::Real`: Standard deviation of the Gaussian blur.

# Returns
- `LF{T}`: Blurred LatticeField.
"""
function lfBlur(f::LF{T,S,N},r::Real) where {T,S<:Complex,N}
    ker = lfGaussian(Intensity,f.L,r).data
    data = isft(sft(ker) .* sft(f.data))
    return LF{T}(data,f.L,f.flambda)
end

function lfBlur(f::LF{T,S,N},r::Real) where {T,S<:Real,N}
    ker = lfGaussian(Intensity,f.L,r).data
    data = abs.(isft(sft(ker) .* sft(f.data)))
    return LF{T}(data,f.L,f.flambda)
end

################################## Metaprogramming for auto-generating methods #################################
#=
The macro @addTemplateMethods(funcExpr) takes a function expression defining a template function (e.g. lfGaussian)
and generates two methods for it: one that accepts a Lattice and FieldVal type, and another that accepts an LF object.
The generated methods handle common tasks such as centering the lattice and formatting the output, allowing the
user-defined function to focus solely on the specific computation.  

In defining the function to be passed to the macro, the following rules apply:
- The function name should be the template function name (e.g. lfGaussian).
- The argument(s) should be the peculiar arguments for the function (e.g. radius for lfGaussian). Both positional and keyword arguments are supported.
- The function body should compute the array `p` representing the field values, using the centered lattice `Lc`.
- The function should not include return statements; the macro will append a standard return statement.
=#

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

################################## LF emojies #################################

# The functions below are for generating fun shapes: Heart, smiley face, and pointing hand emojies. 

#------- Heart -------#
function heartQ(x::Real,y::Real,w::Real,t::Real,b::Real)
    r = 1 - cbrt( (x/w)^2 )
    if r < 0 
        return 0.0
    end
    c = sqrt(r)
    yt = t*(2*c - c^2 - c^3)
    yb = b*(-2*c - c^2 + c^3)
    return yb <= y <= yt ? 1.0 : 0.0
end
function heartQ(x::Real,y::Real,scale::Real)
    w = scale
    t = scale * 1.2/sqrt(2)
    b = scale/sqrt(2)
    return heartQ(x,y,w,t,b)
end

@addTemplateMethods function lfHeart(scale; flip=false)
    p = [heartQ(x,y,scale) for x in Lc[1], y in Lc[2]]
    if flip
        p = reverse(transpose(p),dims=1)
    end
end

#------- Smiley -------#

function smileQ(x::Real,y::Real,hr::Real,mr::Real,ma::Real,mt::Real,erx::Real,ery::Real,ex::Real,ey::Real)
    if x^2+y^2>hr^2   # Outside head
        return 0.0
    end
    if y < 0 && (mr+mt)^2 > x^2+y^2 > (mr-mt)^2 && abs(x/y) < tan(ma)    # In lips
        return 0.1
    end
    if (abs(x) - ex)^2/erx^2 + (y-ey)^2/ery^2 < 1    # In eyes
        return 0.1
    end
    return 1.0
end

function smileQ(x::Real,y::Real,scale::Real)
    hr = scale
    mr = scale * 0.6
    ma = 3*pi/8
    mt = scale * 0.05
    erx = scale * 0.12
    ery = scale * 0.25
    ex = scale * 0.3
    ey = scale * 0.3
    return smileQ(x,y,hr,mr,ma,mt,erx,ery,ex,ey)
end

@addTemplateMethods function lfSmile(scale; flip=false)
    p = [smileQ(x,y,scale) for x in Lc[1], y in Lc[2]]
    if flip
        p = reverse(transpose(p),dims=1)
    end
end

#------- Pointer -------#

function pointerOutlineQ(x::Real,y::Real,bt::Real,bh::Real,
        fr::Real,l1::Real,l2::Real,l3::Real,l4::Real,
        tr::Real,l0::Real,
        hl0::Real,hl1::Real)
    R = (fr*8+tr)/2
    # base of hand
    if x<=0 && (R-bt)^2 <= x^2+y^2 <= (R+bt)^2
        return bh
    end
    # outside of hand
    if (-R-bt)<=y<=(-R+bt) && 0<=x<=hl1
        return bh
    end
    # finger tips
    tcx = hl1 .+ [l1,l2,l3,l4]
    tcy = (-R + fr) .+ 2*fr*[3:-1:0;]
    if any( (fr-bt)^2 <= (x-tcx[j])^2 + (y-tcy[j])^2 <= (fr+bt)^2 && x>tcx[j] for j=1:4)
        return bh
    end

    # Between fingers
    if any( tcy[j]-fr-bt<=y<=tcy[j]-fr+bt && hl1 <= x <= tcx[j] for j=1:4 )
        return bh
    end
    if any( (y-(tcy[j]-fr))^2 + (x-hl1)^2 <= bt^2 for j=1:3 )
        return bh
    end
    
    # Top of hand
    if tcy[1] + fr - bt <= y <= tcy[1] + fr + bt && hl0 <= x <= tcx[1]
        return bh
    end
    if (x-hl0)^2 + (y- (tcy[1]+fr))^2 <= bt^2
        return bh
    end

    # Thumb
    if x>=(hl0+l0) && y >= tcy[1]+fr && (tr-bt)^2 <= (x - (hl0+l0))^2 + (y-(tcy[1]+fr))^2 <= (tr+bt)^2
        return bh
    end
    if hl0+l0>=x>=0 && (R-bt)<=y<=(R+bt)
        return bh
    end
    return 0.0
end

function pointerOutlineQ(x::Real,y::Real,scale::Real)
    bt = scale * 0.025
    bh = 0.5
    fr = scale * 0.1
    l1 = scale * 0.7
    l2 = scale * 0.25
    l3 = scale * 0.2
    l4 = scale * 0.15
    tr = scale * 0.15
    l0 = scale * 0.1
    hl0 = scale * 0.2
    hl1 = scale * 0.5
    return pointerOutlineQ(x,y,bt,bh,fr,l1,l2,l3,l4,tr,l0,hl0,hl1)
end

function pointerFillQ(x::Real,y::Real,
        fr::Real,l1::Real,l2::Real,l3::Real,l4::Real,
        tr::Real,l0::Real,
        hl0::Real,hl1::Real)
    R = (fr*8+tr)/2
    tcx = hl1 .+ [l1,l2,l3,l4]
    tcy = (-R + fr) .+ 2*fr*[3:-1:0;]
    
    # base of hand
    if x<=0 && x^2+y^2 <= R^2
        return 1.0
    end

    # palm
    if -R<=y<=tcy[1]+fr && 0<=x<=hl1
        return 1.0
    end

    # fingers
    if any( tcy[j]-fr <= y <= tcy[j]+fr && hl1<=x<=tcx[j] for j=1:4)
        return 1.0
    end
    if any( x>tcx[j] && (x-tcx[j])^2 + (y-tcy[j])^2 <= fr^2 for j=1:4)
        return 1.0
    end

    # thumb
    if 0<=x<=hl0+l0 && tcy[1]+fr<=y<=R
        return 1.0
    end
    if x>=hl0+l0 && y>=tcy[1]+fr && (x-(hl0+l0))^2 + (y-(tcy[1]+fr))^2 < tr^2
        return 1.0
    end
    return 0.0
end

function pointerFillQ(x::Real,y::Real,scale::Real)
    fr = scale * 0.1
    l1 = scale * 0.7
    l2 = scale * 0.25
    l3 = scale * 0.2
    l4 = scale * 0.15
    tr = scale * 0.15
    l0 = scale * 0.1
    hl0 = scale * 0.2
    hl1 = scale * 0.5
    return pointerFillQ(x,y,fr,l1,l2,l3,l4,tr,l0,hl0,hl1)
end

@addTemplateMethods function lfPointer(scale; flip=false)
    dataFill = [pointerFillQ(x,y,scale) for x in Lc[1], y in Lc[2]]
    dataBorder = [pointerOutlineQ(x,y,scale) for x in Lc[1], y in Lc[2]]
    p = [iszero(dataBorder[J]) ? dataFill[J] : dataBorder[J] for J in CartesianIndices(size(dataFill))]
    if flip
        p = reverse(transpose(p),dims=1)
    end
end