module LatticeCore
using FFTW: fftshift, ifftshift
export Lattice, natlat, padout, sft, isft, elq
export RealPhase, Generic, Phase, FieldVal, ComplexPhase

"""
    Lattice{N}

    A type alias for representing an N-dimensional lattice.

    Each element in the `NTuple` is an `AbstractRange`, representing a range in each dimension of the lattice.

    # Type Parameters
    - `N`: The number of dimensions in the lattice.

    """
const Lattice{N} = NTuple{N,AbstractRange} where {N}

"""
   elq(x::Lattice{N}, y::Lattice{N}) where {N} 
   
    Checks if two lattices are equal, throwing a `DomainError` if they are not.
"""
elq(x::Lattice{N}, y::Lattice{N}) where {N} = (all(isapprox.(x, y)) || throw(DomainError((x, y), "Unequal lattices.")); return nothing)



"""
    FieldVal

Abstract type representing the general category of field values in spatial light modulator (SLM) fields. 
It serves as a base type for various specific types of field values such as phase, intensity, and amplitude values. 
Subtypes of `FieldVal` are used to define and work with different properties of SLM fields in a structured manner.

# Examples
- `Phase <: FieldVal`
- `Intensity <: FieldVal`
- `Amplitude <: FieldVal`
"""
abstract type FieldVal end

"""
    Generic <: FieldVal

A subtype of `FieldVal`, representing a generic field value. This type can be used as a placeholder or 
for cases where the specific nature of the field value is not critical to the operation or function being performed. 
It allows for flexibility and can be replaced or specified further in more detailed implementations.

# Usage
- Useful in functions or algorithms where the specific type of field value is not a primary concern.
"""
abstract type Generic <: FieldVal end

"""
    Phase <: FieldVal

Abstract type representing phase values in SLM fields. `Phase` is a subtype of `FieldVal` and acts as 
a parent type for different kinds of phase representations like real and complex phases. It is central 
to operations and calculations involving phase modulation and manipulation in optical systems.

# Subtypes
- `RealPhase <: Phase`
- `ComplexPhase <: Phase`
"""
abstract type Phase <: FieldVal end

"""
    RealPhase <: Phase

A subtype of `Phase`, representing real phase values in SLM fields. `RealPhase` is used for scenarios 
where the phase values are real numbers. This is commonly encountered in various optical applications 
where phase modulation does not involve complex numbers.

# See Also
- `ComplexPhase <: Phase` for complex phase values.
"""
abstract type RealPhase <: Phase end



"""
    const UPhase = RealPhase

Type alias for `RealPhase`, representing a specific kind of real phase value in SLM fields. `UPhase` is typically 
used to denote a standard real phase, possibly indicating a certain format or representation used in the system.
It inherits all properties and characteristics of `RealPhase`.

# Usage
- Use `UPhase` in place of `RealPhase` for clearer code semantics or when referring to a standard format of real phase.
"""
const UPhase = RealPhase

"""
    const UnwrappedPhase = RealPhase

Type alias for `RealPhase`, specifically used to represent real phase values that have been unwrapped. 
`UnwrappedPhase` indicates that the phase values are processed to remove discontinuities, a common procedure 
in phase analysis. This alias helps in explicitly indicating the nature of the phase data being handled.

# Usage
- Use `UnwrappedPhase` to signify that the phase values have undergone unwrapping.
"""
const UnwrappedPhase = RealPhase

"""
    ComplexPhase <: Phase

Abstract type representing complex phase values in SLM fields. `ComplexPhase` is a subtype of `Phase` and is used 
in contexts where phase values are complex numbers. This type is essential for operations that involve complex 
phase modulation, such as in certain advanced optical trapping setups.

# See Also
- `RealPhase <: Phase` for real phase values.
"""
abstract type ComplexPhase <: Phase end

"""
    const S1Phase = ComplexPhase

Type alias for `ComplexPhase`, representing a specific format or representation of complex phase values in SLM fields. 
`S1Phase` might indicate a certain standard or convention in complex phase handling. It is used to provide clarity 
and specificity in code where complex phases are used in a particular context or format.

# Usage
- Use `S1Phase` to specify the use of a particular format or standard of complex phase values.
"""
const S1Phase = ComplexPhase





""" 
    natlat(n::Int) 
    Create a "natural" lattice for a given integer n.
    This funnction generates a 1D lattice with n points, ranging from -floor(n/2) to floor((n-1)/2),
    and scales it by 1/sqrt(n) to normalize.
    The natural lattice is often used in FFT problems to represent the discrete points in real space."""
function natlat(n::Int)
    return (-floor(n / 2):floor((n - 1) / 2)) ./ sqrt(n)
end

"""
    natlat(s::NTuple{N,Integer}) 
    Overload of natlat for an N-dimensional lattice specified by a tuple of integers."""
function natlat(s::NTuple{N,Integer}) where {N}
    return ((natlat(n) for n in s)...,)
end

"""
    padout(L::AbstractRange, p::Tuple{Int,Int}) 
    Pads an AbstractRange by a certain number of points before the first element and after the last element.
    The padding amount is specified by the tuple `p`, where `p[1]` is the padding before and `p[2]` is the padding after.
    It returns a new range with the specified padding added, scaled by the step of the original range,
    and shifted to start at the correct value. Several overloads are provided for convenience."""
function padout(L::AbstractRange, p::Tuple{Int,Int})
    return (0:sum(p)+length(L)-1) .* step(L) .+ (L[1] - step(L) * p[1])
end

#then some syntactic sugar for padout for abstract ranges
padout(L::AbstractRange, p::Int) = padout(L, (p, p))
padout(L::Lattice{N}, p::NTuple{M,Int}) where {N,M} = M == 2N ? ((padout(L[i], (p[2i-1], p[2i])) for i = 1:N)...,) : error("Bad tuple length.")
padout(L::Lattice{N}, p::NTuple{N,Int}) where {N} = ((padout(L[i], p[i]) for i = 1:N)...,)
padout(L::Lattice{N}, p::Int) where {N} = ((padout(L[i], p) for i = 1:N)...,)

"""Pads an N-dimensional array `A` with a certain padding specified by the tuple `p` and uses `filler` as the padding value.
    The `filler` defaults to zero of the same type as the array. Several overloads are provided for convenience.
    """
function padout(A::Array{T,N}, p::NTuple{M,Int}, filler=zero(T)) where {T,N,M}
    M == 2N || error("Bad tuple length.")
    output = fill(filler, size(A) .+ ((p[2i-1] + p[2i] for i = 1:N)...,))
    output[((p[2i-1]+1:size(output, i)-p[2i]) for i = 1:N)...] .= A
    return output
end
#then some syntactic sugar for padout for arrays
padout(A::Array{T,N}, p::NTuple{N,Int}, filler=zero(T)) where {T,N} = padout(A, ((p[ceil(Int, i / 2)] for i = 1:2N)...,), filler)
padout(A::Array{T,N}, p::Int, filler=zero(T)) where {T,N} = padout(A, ((p for i = 1:N)...,), filler)


"""
    sft(v)

    Performs a shifted Fourier transform on the input vector `v`. This function first applies an inverse shift (using `ifftshift`), then performs a Fourier transform (`fft`), and finally applies a forward shift (`fftshift`).
    # Arguments
    - `v`: An array to be transformed.
    # Returns
    - An array containing the shifted Fourier transform of `v`.
    """
sft(v) = fftshift(fft(ifftshift(v)))


"""
    isft(v)
    Performs an inverse shifted Fourier transform on the input vector `v`. This function first applies an inverse shift (using `ifftshift`), then performs an inverse Fourier transform (`ifft`), and finally applies a forward shift (`fftshift`).
    # Arguments
    - `v`: An array to be transformed.
    # Returns
    - An array containing the inverse shifted Fourier transform of `v`.
    """
isft(v) = fftshift(ifft(ifftshift(v)))




end # module Core