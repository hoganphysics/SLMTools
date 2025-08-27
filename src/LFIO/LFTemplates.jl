#=
Template functions for generating common LatticeField distributions.  These are implemented with metaprogramming to reduce boilerplate.

Requires:
    using FreeTypeAbstraction: findfont, renderstring!
=#

export lfParabola, lfGaussian, lfRing, lfCap, ftaText, lfText, lfRect, lfRand

################################## Function docstrings #################################

"""
    lfParabola(quad::Real, lin::NTuple{N,Real}=Tuple(0.0 for i=1:N))

    LF template function generating a parabolic LatticeField with specified quadratic and linear coefficients.  All 
    LF template functions have a method which accepts a Lattice and LF type (e.g. Intensity, Modulus) and another
    method which acccepts an LF from which the Lattice and LF type are inferred. In either case, the first argument
    or arguments are these template objects, and subsequent arguments are peculiar to the function.  See examples 
    below for how to use the template arguments.  In this case, the peculiar arguments are as follows: 

    # Peculiar arguments
    - `quad::Real`: Second derivative of the parabola.
    - `lin::NTuple{N,Real}`: First derivative of the parabola in each dimension.

    # Returns
    - `LatticeField` with a parabolic distribution.

    # Notes
    - The output is calculated as `quad * r2(L)/2 + ldot(lin,L)`, where `r2(L)` is the squared radius
      in the lattice coordinates and `ldot` is the lattice dot product.

    # Examples
    - f1 = lfParabola(Intensity,natlat(128,128),1.5)
    - f2 = lfParabola(f1,2.1,(0.1,0.2))
    """
function lfParabola()
    println("This is a placeholder method for documenting lfParabola. See docstring for functional methods.")
end

"""
    lfGaussian(radius::Real, norm::Real=1.0)

    LF template function generating a Gaussian LatticeField with specified radius and normalization.  All 
    LF template functions have a method which accepts a Lattice and LF type (e.g. Intensity, Modulus) and another
    method which acccepts an LF from which the Lattice and LF type are inferred. In either case, the first argument
    or arguments are these template objects, and subsequent arguments are peculiar to the function.  See examples 
    below for how to use the template arguments.  In this case, the peculiar arguments are as follows: 

    # Peculiar arguments
    - `radius::Real`: Standard deviation of the Gaussian.
    - `norm::Real`: Normalization factor for the Gaussian.

    # Returns
    - `LatticeField` with a Gaussian distribution.

    # Notes
    - The output is calculated as `exp(-r2(L)/(2*radius^2))`, where `r2(L)` is the squared radius
      in the lattice coordinates.

    # Examples
    - f1 = lfGaussian(Intensity,natlat(128,128),1.5)
    - f2 = lfGaussian(f1,2.1)
    """
function lfGaussian()
    println("This is a placeholder method for documenting lfGaussian. See docstring for functional methods.")
end

"""
    lfRing(radius::Number, width::Number)

    LF template function generating a ring-shaped LatticeField with specified radius and width.  All 
    LF template functions have a method which accepts a Lattice and LF type (e.g. Intensity, Modulus) and another
    method which acccepts an LF from which the Lattice and LF type are inferred. In either case, the first argument
    or arguments are these template objects, and subsequent arguments are peculiar to the function.  See examples 
    below for how to use the template arguments.  In this case, the peculiar arguments are as follows: 

    # Peculiar arguments
    - `radius::Number`: Radius of the ring.
    - `width::Number`: Width of the ring.

    # Returns
    - `LatticeField` with a ring-shaped distribution.

    # Notes
    - The output is calculated as `exp(-(sqrt(r2(L)) - radius)^2 / (2*width^2))`, where `r2(L)` is the squared radius
      in the lattice coordinates.

    # Examples
    - f1 = lfRing(Intensity,natlat(128,128),2.0,0.5)
    - f2 = lfRing(f1,2.0,0.5)
    """
function lfRing()
    println("This is a placeholder method for documenting lfRing. See docstring for functional methods.")
end

    """
    lfCap(curvature::Real, height::Real)

    LF template function generating a truncated parabolic LatticeField with specified curvature and height.  All 
    LF template functions have a method which accepts a Lattice and LF type (e.g. Intensity, Modulus) and another
    method which acccepts an LF from which the Lattice and LF type are inferred. In either case, the first argument
    or arguments are these template objects, and subsequent arguments are peculiar to the function.  See examples 
    below for how to use the template arguments.  In this case, the peculiar arguments are as follows: 

    # Peculiar arguments
    - `curvature::Real`: Second derivative of the parabola.
    - `height::Real`: Height of the cap.

    # Returns
    - `LatticeField` with a parabola cap distribution.

    # Notes
    - The output is calculated as `ramp(curvature * r2(L)/2)`, where `r2(L)` is the squared radius
      in the lattice coordinates.

    # Examples
    - f1 = lfCap(Intensity,natlat(128,128),1.5,2.0)
    - f2 = lfCap(f1,2.1,2.0)
    """
    function lfCap()
        println("This is a placeholder method for documenting lfCap. See docstring for functional methods.")
    end

    """
    lfText(str::String)

    LF template function generating a text LatticeField with specified string.  All 
    LF template functions have a method which accepts a Lattice and LF type (e.g. Intensity, Modulus) and another
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

    # Returns
    - `LatticeField` with a text distribution.

    # Notes
    - The output is generated using the `ftaText` function, which renders the text string into a float array.

    # Examples
    - f1 = lfText(Intensity,natlat(128,128),"C")
    - f2 = lfText(f1,"C")
    """
    function lfText()
        println("This is a placeholder method for documenting lfText. See docstring for functional methods.")
    end

    """
    lfRect()

    LF template function generating a rectangular shaped (i.e. boxcar function) LatticeField.  All
    LF template functions have a method which accepts a Lattice and LF type (e.g. Intensity, Modulus) and another
    method which acccepts an LF from which the Lattice and LF type are inferred. In either case, the first argument
    or arguments are these template objects, and subsequent arguments are peculiar to the function.  See examples 
    below for how to use the template arguments.  In this case, the peculiar arguments are as follows: 
    
    # Peculiar arguments
    - sides::NTuple{N,Real}: Lengths of the sides of the rectangle in each dimension.
    - height::Real=1.0: Height of the boxcar function.

    # Returns
    - `LatticeField` with a rectangular shape.

    # Examples
    - f1 = lfRect(Intensity,natlat(128,128),(1.5,2.3))
    - f2 = lfRect(f1,(1.5,1.3))
    """
    function lfRect()
        println("This is a placeholder method for documenting lfRect. See docstring for functional methods.")
    end

    """
    lfRand()

    LF template function generating a random LatticeField.  All
    LF template functions have a method which accepts a Lattice and LF type (e.g. Intensity, Modulus) and another
    method which acccepts an LF from which the Lattice and LF type are inferred. In either case, the first argument
    or arguments are these template objects, and subsequent arguments are peculiar to the function.  See examples 
    below for how to use the template arguments.  In this case, there are no peculiar arguments, save an optional 
    keyword argument to specify the data type of the field values.

    # Peculiar arguments
    - R::DataType=Float64: Data type of the random values.
    
    # Returns
    - `LatticeField` with random values.

    # Examples
    - f1 = lfRand(Intensity,natlat(128,128))
    - f2 = lfRand(f1)
    """
    function lfRand()
        println("This is a placeholder method for documenting lfRand. See docstring for functional methods.")
    end

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
