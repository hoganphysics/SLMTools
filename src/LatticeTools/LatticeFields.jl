module LatticeFields

export Lattice, elq, RealPhase, Generic, Phase, FieldVal, ComplexPhase, UPhase, UnwrappedPhase, S1Phase, Intensity, Amplitude, Modulus, RealAmplitude, RealAmp, ComplexAmplitude, ComplexAmp, LatticeField, LF, normalizeLF
export subfield, wrap, square, sublattice

#region ---------------------------Abstract Types and aliases--------------------------------------------
"""
    Lattice{N}

    A type alias for representing an N-dimensional lattice.

    Each element in the `NTuple` is an `AbstractRange`, representing a range in each dimension of the lattice.

    # Type Parameters
    - `N`: The number of dimensions in the lattice.

    """
const Lattice{N} = NTuple{N,AbstractRange} where {N}

"""
    FieldVal

    Abstract supertype for all types used in a LatticeField. 
    It serves as a parent type for various other abstract types like Phase, Intensity, Modulus, etc. 
    Subtypes of `FieldVal` are used to provide customized dispatch depending on what a LatticeField represents.

    # Examples
    - `Phase <: FieldVal`
    - `Intensity <: FieldVal`
    - `Amplitude <: FieldVal`
    """
abstract type FieldVal end

"""
    Generic <: FieldVal

    A subtype of `FieldVal`, representing a generic field value. Placeholder type most for testing.
    """
abstract type Generic <: FieldVal end

"""
    Phase <: FieldVal

    Abstract type representing phase values in SLM fields. `Phase` is a subtype of `FieldVal` and acts as 
    a parent type for RealPhase and ComplexPhase. 

    # Subtypes
    - `RealPhase <: Phase`
    - `ComplexPhase <: Phase`
    """
abstract type Phase <: FieldVal end

"""
    RealPhase <: Phase

    A subtype of `Phase`, representing an "unwrapped" real-valued phase.

    # See Also
    - `ComplexPhase <: Phase` for "wrapped" complex-valued phases.
    """
abstract type RealPhase <: Phase end

"""
    const UPhase = RealPhase

    Type alias for `RealPhase`, short for "unwrapped phase".
    """
const UPhase = RealPhase

"""
    const UnwrappedPhase = RealPhase

    Type alias for `RealPhase`.
    """
const UnwrappedPhase = RealPhase

"""
    ComplexPhase <: Phase

    Abstract type representing "wrapped" complex-valued phases. 

    # See Also
    - `RealPhase <: Phase` for "unwrapped" real-valued phases.
    """
abstract type ComplexPhase <: Phase end

"""
    const S1Phase = ComplexPhase

    Type alias for `ComplexPhase`, "S1" referring to the complex unit circle where these phases live.
    """
const S1Phase = ComplexPhase

"""
    Intensity <: FieldVal

    Abstract type representing intensity values in SLM fields.

    # Usage
    - Use `Intensity` as a base type for defining various forms of light intensity data.
    """
abstract type Intensity <: FieldVal end

"""
    Amplitude <: FieldVal

    Abstract type representing amplitude values. Parent type to Modulus and ComplexAmplitude, the difference between the two
	being that the latter includes phase information while the former does not. 

    # Subtypes
    - `Modulus <: Amplitude` for non-negative real amplitudes.
    - `ComplexAmplitude <: Amplitude` for complex amplitudes.
    """
abstract type Amplitude <: FieldVal end

"""
    Modulus <: Amplitude

    Abstract type representing non-negative real amplitude values.  Subtype of Amplitude.

    # See Also
    - `ComplexAmplitude <: Amplitude` for complex amplitude representations.
    """
abstract type Modulus <: Amplitude end

"""
    const RealAmplitude = Modulus

    Type alias for `Modulus`.
    """
const RealAmplitude = Modulus

"""
    const RealAmp = Modulus

    Type alias for `Modulus`.
    """
const RealAmp = Modulus

"""
    ComplexAmplitude <: Amplitude

    Abstract type representing complex amplitude values.  Subtype of `Amplitude`.  Incorporates both phase and magnitude information.  
	Must be used in conjunction with a Complex numerical data type. 

    # See Also
    - `Modulus <: Amplitude` for non-negative real amplitude representations.
    """
abstract type ComplexAmplitude <: Amplitude end

"""
        const ComplexAmp = ComplexAmplitude

    Type alias for `ComplexAmplitude`.
    """
const ComplexAmp = ComplexAmplitude
#endregion

#region ---------------------------LatticeField--------------------------------------------

"""
    LatticeField{S<:FieldVal,T,N}

    A struct representing a field defined on a lattice. It encapsulates the field data along with its lattice structure and an optional scaling factor for the wavelength*focal length product.

    # Fields
    - `data::AbstractArray{T,N}`: The field data, stored as an N-dimensional array.
    - `L::Lattice{N}`: The lattice structure, defined as an N-tuple of ranges.
    - `flambda::Real`: A scaling factor for the wavelength*focal length product, with a default value of 1.0.

    # Parameters
    - `S`: The subtype of `FieldVal` representing the type of field value (e.g., RealPhase, Intensity).
    - `T`: The element type of the field data array.
    - `N`: The dimensionality of the lattice and data array.

    # Constructor
    The constructor checks that the size of the `data` array matches the size of the lattice `L`. If they do not match, 
    a `DimensionMismatch` error is thrown.

    # Examples
    data = exp.(2*pi*im*rand(Float64, 10, 10))
    lattice = (1:10, 1:10)
    lf = LatticeField{ComplexPhase}(data, lattice)
    This will create a LatticeField instance with complex phase values on a 10x10 lattice.

    """
struct LatticeField{S<:FieldVal,T,N}
    data::AbstractArray{T,N}
    L::Lattice{N}
    flambda::Real
    function LatticeField{S,T,N}(data::AbstractArray{T,N}, L::Lattice{N}, flambda::Real=1.0) where {S<:FieldVal,T,N}
        if size(data) !== length.(L)
            throw(DimensionMismatch("Field data size does not match lattice size."))
        else
            return new(data, L, flambda)
        end
    end
end

"""
    LatticeField{S}(array::AbstractArray{T,N}, L::Lattice{N}, flambda::Real=1.0) where {S<:FieldVal,T,N}

    Convenience constructor for `LatticeField`. It allows for creating an instance of `LatticeField` by 
    specifying the type of field value `S`, and inferring the data type `T` and dimension `N` from the provided array. 
    The `flambda` parameter is an optional scaling factor for the wavelength, with a default value of 1.0.
	For Intensity fields, the input array is ramped to ensure it is non-negative.  This avoid errors when converting to Modulus. 

    # Parameters
    - `array`: The field data as an N-dimensional array.
    - `L`: The lattice structure, represented as an N-tuple of ranges.
    - `flambda`: Wavelength*focal length scale factor (optional).  Default value == 1.0.
    - `S`: A subtype of `FieldVal` representing the field value type.
    - `T`: Inferred from the element type of `array`.
    - `N`: Inferred from the dimensionality of `array`.

    This constructor is particularly useful when the types `T` and `N` can be automatically inferred from the input array. 
    """
LatticeField{S}(array::AbstractArray{T,N}, L::Lattice{N}, flambda::Real=1.0) where {S<:FieldVal,T,N} = LatticeField{S,T,N}(array, L, flambda)
ramp(x::T) where {T<:Number} = (x < 0 ? zero(T) : x)
LatticeField{Intensity}(array::AbstractArray{T,N},L::Lattice{N},flambda::Real=1.0) where {T,N} = LatticeField{Intensity,T,N}(ramp.(array),L,flambda)

"""
    const LF = LatticeField

    Type alias for the `LatticeField` struct.
    """
const LF = LatticeField

"""
   elq(x::Lattice{N}, y::Lattice{N}) where {N} 
   
    Checks if two lattices are equal, throwing a `DomainError` if they are not.
    """
elq(x::Lattice{N}, y::Lattice{N}) where {N} = (all(isapprox.(x, y)) || throw(DomainError((x, y), "Unequal lattices.")); return nothing)

# Equal lattice query
elq(x::LF, y::LF) = (all(isapprox.(x.L, y.L)) && isapprox(x.flambda,y.flambda) || throw(DomainError((x.L, y.L), "Unequal lattices.")); return nothing)


# Minimal methods to make a LatticeField function similarly to an array.
Base.size(f::LatticeField) = size(f.data)
Base.getindex(f::LatticeField, i::Int) = getindex(f.data, i)
Base.setindex!(f::LatticeField, v, i::Int) = setindex!(f.data, v, i)
Base.setindex!(f::LatticeField{S,T,N}, v, indices::Vararg{Int,N}) where {S,T,N} = setindex!(f.data, v, indices...)
Base.IndexStyle(::Type{<:LatticeField}) = IndexLinear()
Base.axes(f::LatticeField, d) = axes(f.data, d)
Base.show(io::IO, f::T) where {T<:LF} = (println(T); println("Lattice: ", f.L); display(f.data))



#= 
# The following functions allow for creation of sublattice fields by subscripting. 
=#
function Base.getindex(f::LatticeField{S,T,N}, x::I...) where {S,T,N,I<:Integer}
    # Subscripting f by integers returns the field value at that location. 
    return getindex(f.data, x...)
end
function Base.getindex(f::LatticeField{S,T,N}, x::AbstractRange{I}...) where {S,T,N,I<:Integer}
    # Subscripting f by ranges returns a subfield on the sublattice. 
    return LatticeField{S,T,N}(getindex(f.data, x...), sublattice(f.L, x...), f.flambda)
end
function Base.getindex(f::LatticeField{S,T,N}, x::Union{I,AbstractRange{I}}...) where {S,T,N,I<:Integer}
    # Subscripting f by integers and ranges returns a subfield on the slice fixing dimensions with integer indices. 
    # Note that this is not a type safe method, so it should not be used in performance demanding situations. 
    z = filter(i -> !(x[i] isa Integer), 1:length(x))
    return LatticeField{S,T,length(z)}(getindex(f.data, x...), sublattice(f.L[z], x[z]...), f.flambda)
end
function Base.getindex(f::LatticeField{S,T,N}, x::CartesianIndices) where {S,T,N}
    # Subscripting f by CartesianIndices returns the subfield on that region. 
    return LatticeField{S,T,N}(getindex(f.data,x),sublattice(f.L,x),f.flambda)
end

#endregion

#region ----------------------sublattice--------------------------------
"""
    sublattice(L::Lattice{N}, box::CartesianIndices) where {N}

    Extracts a sublattice from a given `Lattice` `L` using the specified `box` of Cartesian indices.

    # Arguments
    - `L`: A lattice from which the sublattice is to be extracted.
    - `box`: Cartesian indices defining the region of the sublattice.

    # Returns
    - A tuple representing the extracted sublattice.
    - Throws an error if the dimensions of the `box` do not match the lattice.
    """
function sublattice(L::Lattice{N}, box::CartesianIndices) where {N}
    length(size(box)) == N || error("box should have same number of dimensions as lattice.") #TODO check if this should be here
    return ((L[j][box[1][j]:box[end][j]] for j = 1:N)...,)
end

"""
    sublattice(L::Lattice{N}, x::Union{I,AbstractRange{I}}...) where {N,I<:Integer}

    Create a sublattice from a given lattice `L` by specifying indices or ranges in each dimension. 

    # Arguments
    - `L::Lattice{N}`: The original lattice from which the sublattice is to be created.
    - `x::Union{I,AbstractRange{I}}...`: A series of integers or ranges specifying the indices in each dimension to create the sublattice.

    # Returns
    - `Tuple`: A tuple representing the sublattice.

    # Errors
    Throws `DimensionMismatch` if the number of indices provided does not match the dimensions of `L`.

    """
function sublattice(L::Lattice{N}, x::Union{I,AbstractRange{I}}...) where {N,I<:Integer}
    length(L) != length(x) && throw(DimensionMismatch("Wrong number of indices while attempting to make sublattice."))
    y = (((i isa Integer) ? (i:i) : i for i in x)...,)
    return ((L[j][y[j]] for j = 1:N)...,)
end
#endregion

#region ---------------------------subfield--------------------------------------------
"""
    subfield(f::LatticeField{S,T,N}, x::I...) where {S,T,N,I<:Integer}

    Extracts the field value at a specific point in the `LatticeField`. This function allows for subscripting 
    a `LatticeField` instance by integer indices to access individual field values.

    # Parameters
    - `f`: An instance of `LatticeField`.
    - `x`: A series of integer indices corresponding to the point in the lattice.

    # Returns
    The field value at the specified point in the lattice.

    This function is typically used for accessing values at specific lattice points.
    """
function subfield(f::LatticeField{S,T,N}, x::I...) where {S,T,N,I<:Integer}
    # Subscripting f by integers returns the field value at that location. 
    return getindex(f.data, x...)
end

"""
    subfield(f::LatticeField{S,T,N}, x::AbstractRange{I}...) where {S,T,N,I<:Integer}

    Creates a new `LatticeField` representing a subfield of the original field, defined on a sublattice. 
    This function is used for subscripting a `LatticeField` by ranges to extract a portion of the field.

    # Parameters
    - `f`: An instance of `LatticeField`.
    - `x`: A series of ranges specifying the portion of the lattice to extract.

    # Returns
    A new `LatticeField` instance representing the specified subfield.

    This is useful for operations that require working with a specific region of the lattice.
    """
function subfield(f::LatticeField{S,T,N}, x::AbstractRange{I}...) where {S,T,N,I<:Integer}
    # Subscripting f by ranges returns a subfield on the sublattice. 
    return LatticeField{S,T,N}(getindex(f.data, x...), sublattice(f.L, x...), f.flambda)
end

"""
    subfield(f::LatticeField{S,T,N}, x::Union{I,AbstractRange{I}}...) where {S,T,N,I<:Integer}

    Extracts a subfield from a `LatticeField`, allowing for a mix of integer indices and ranges. This function 
    enables subscripting by both integers and ranges, providing flexibility in specifying the subfield.

    # Parameters
    - `f`: An instance of `LatticeField`.
    - `x`: A series of integers and/or ranges specifying the dimensions and regions of the lattice to include.

    # Returns
    A new `LatticeField` instance representing the specified subfield.

    # Note
    This method is not type safe and might not be suitable for performance-critical situations. It is intended 
    for more flexible data extraction where performance is not the primary concern.

    This function is especially useful when needing to fix certain dimensions while varying others.
"""
function subfield(f::LatticeField{S,T,N}, x::Union{I,AbstractRange{I}}...) where {S,T,N,I<:Integer}
    # Subscripting f by integers and ranges returns a subfield on the slice fixing dimensions with integer indices. 
    # Note that this is not a type safe method, so it should not be used in performance demanding situations. 
    z = filter(i -> !(x[i] isa Integer), 1:length(x))
    return LatticeField{S,T,length(z)}(getindex(f.data, x...), sublattice(f.L[z], x[z]...), f.flambda)
end
#endregion

#region -------------------- Basic arithmetic for LatticeFields -------------------------------

# Default behavior is to throw an undefined method error.
Base.:(*)(args::LF...) = error("Behavior undefined for this combination of inputs: ", *((string(i) * ", " for i in typeof.(args))...))
Base.:(+)(args::LF...) = error("Behavior undefined for this combination of inputs: ", *((string(i) * ", " for i in typeof.(args))...))
Base.:(-)(args::LF...) = error("Behavior undefined for this combination of inputs.", *((string(i) * ", " for i in typeof.(args))...))
Base.:(/)(args::LF...) = error("Behavior undefined for this combination of inputs.", *((string(i) * ", " for i in typeof.(args))...))
Base.sqrt(args::LF...) = error("Behavior undefined for this combination of inputs.", *((string(i) * ", " for i in typeof.(args))...))
square(args::LF...) = error("Behavior undefined for this combination of inputs.", *((string(i) * ", " for i in typeof.(args))...))
Base.abs(args::LF...) = error("Behavior undefined for this combination of inputs.", *((string(i) * ", " for i in typeof.(args))...))

# Intensity to Modulus
Base.sqrt(f::LF{Intensity}) = LF{Modulus}(sqrt.(f.data), f.L, f.flambda)

# Amplitude times phase
Base.:(*)(x::LF{<:Amplitude}, y::LF{ComplexPhase}) = (elq(x, y); LF{ComplexAmp}(x.data .* y.data, x.L, x.flambda))
Base.:(*)(x::LF{<:Amplitude}, y::LF{RealPhase}) = (elq(x, y); LF{ComplexAmp}(x.data .* exp.(2pi * im * y.data), x.L, x.flambda))
Base.:(*)(y::LF{ComplexPhase}, x::LF{<:Amplitude}) = (elq(x, y); LF{ComplexAmp}(x.data .* y.data, x.L, x.flambda))
Base.:(*)(y::LF{RealPhase}, x::LF{<:Amplitude}) = (elq(x, y); LF{ComplexAmp}(x.data .* exp.(2pi * im * y.data), x.L, x.flambda))

# Combining phases
Base.:(*)(x::LF{ComplexPhase}, y::LF{ComplexPhase}) = (elq(x, y); LF{ComplexPhase}(x.data .* y.data, x.L, x.flambda))
Base.:(*)(x::LF{RealPhase}, y::LF{ComplexPhase}) = (elq(x, y); LF{ComplexPhase}(exp.(2pi * im * x.data) .* y.data, x.L, x.flambda))
Base.:(*)(x::LF{ComplexPhase}, y::LF{RealPhase}) = (elq(x, y); LF{ComplexPhase}(x.data .* exp.(2pi * im * y.data), x.L, x.flambda))

Base.:(+)(x::LF{RealPhase}, y::LF{RealPhase}) = (elq(x, y); LF{RealPhase}(x.data .+ y.data, x.L, x.flambda))

# ComplexAmplitude to Modulus
Base.abs(x::LF{ComplexAmp}) = LF{Modulus}(abs.(x.data), x.L, x.flambda)

# Phase conjugation
Base.conj(x::LF{RealPhase}) = LF{RealPhase}(-x.data, x.L, x.flambda)
Base.conj(x::LF{ComplexPhase}) = LF{ComplexPhase}(conj.(x.data), x.L, x.flambda)

# Scalar operations
Base.:(*)(x::T, f::LF{S,T,N}) where {S,T,N} = LF{S,T,N}(f.data .* x, f.L, f.flambda)
Base.:(*)(f::LF{S,T,N}, x::T) where {S,T,N} = LF{S,T,N}(f.data .* x, f.L, f.flambda)
Base.:(+)(x::T, f::LF{S,T,N}) where {S,T,N} = LF{S,T,N}(f.data .+ x, f.L, f.flambda)
Base.:(+)(f::LF{S,T,N}, x::T) where {S,T,N} = LF{S,T,N}(f.data .+ x, f.L, f.flambda)
Base.:(/)(f::LF{S,T,N}, x::T) where {S,T,N} = LF{S,T,N}(f.data ./ x, f.L, f.flambda)

#endregion

#region --------------------------- wrap, square, normalizeLF --------------------------------------------
"""
    square(x::LF{<:Amplitude})

    Transforms a `LatticeField` containing amplitude values (`Amplitude` or its subtypes like `Modulus` or `ComplexAmplitude`) 
    into a `LatticeField` of intensity values. The transformation is achieved by squaring the absolute value of each element 
    in the `data` field of the `LatticeField`.

    # Parameters
    - `x`: An instance of `LatticeField` with a subtype of `Amplitude`.

    # Returns
    A new `LatticeField` instance with the `Intensity` type, where each data point is the square of the absolute value 
    of the corresponding point in the input `LatticeField`.
    """
square(x::LF{<:Amplitude}) = LF{Intensity}(abs.(x.data) .^ 2, x.L, x.flambda)

"""
    wrap(x::LF{RealPhase})

    Transforms a `LatticeField` containing real phase values into a `LatticeField` of complex phase values. 
    The transformation involves converting each real phase value to a complex phase value using Euler's formula.

    # Parameters
    - `x`: An instance of `LatticeField` with `RealPhase`.

    # Returns
    A new `LatticeField` instance with `ComplexPhase`, where each data point is the result of `exp(2πi * data)` 
    of the corresponding point in the input `LatticeField`.
    """
wrap(x::LF{RealPhase}) = LF{ComplexPhase}(exp.(2pi * im * x.data), x.L, x.flambda)

"""
    wrap(x::LF{ComplexPhase})

    Returns a copy of the input `LatticeField` if it already contains complex phase values. This version of the `wrap` 
    function is essentially a no-operation for `LatticeField` instances that are already of the `ComplexPhase` type.

    # Parameters
    - `x`: An instance of `LatticeField` with `ComplexPhase`.

    # Returns
    A copy of the input `LatticeField`.

    This function ensures compatibility and ease of use in code where the exact phase type of a `LatticeField` might not 
    be known in advance.
    """
wrap(x::LF{ComplexPhase}) = LF{ComplexPhase}(x.data,x.L,x.flambda)


"""
    normalizeLF(f::LF{Intensity,T,N})

    Normalizes an LF{Intensity} to have sum 1, i.e. to be a probability distribution. 

    # Parameters
    - `f`: An intensity `LatticeField`.

    # Returns
    A normalized intensity `LatticeField`.
    """
function normalizeLF(f::LF{Intensity,T,N}) where {T<:Real,N}
    # Normalizes an LF{Intensity} to have sum 1, i.e. to be a probability distribution. 
    return LF{Intensity,T,N}(abs.(f.data ./ sum(f.data)), f.L,f.flambda)
end
"""
    normalizeLF(f::LF{<:Amplitude,T,N})

    Normalizes an LF{<:Amplitude} so that the corresponding intensity has sum 1, i.e. the intensity is a probability distribution. 

    # Parameters
    - `f`: An intensity `LatticeField`.

    # Returns
    A normalized intensity `LatticeField`.
    """
function normalizeLF(f::LF{S,T,N}) where {S<:Amplitude,T<:Real,N}
    # Normalizes an LF{<:Amplitude} so that the corresponding intensity has sum 1, i.e. the intensity is a probability distribution. 
    return LF{S,T,N}(f.data ./ sqrt(sum(abs.(f.data).^2)), f.L,f.flambda)
end

#endregion



end # module LatticeFields