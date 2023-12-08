module LatticeCore
include("LatticeFields.jl")
using .LatticeFields
export Lattice, elq, RealPhase, Generic, Phase, FieldVal, ComplexPhase, UPhase, UnwrappedPhase, S1Phase
export Intensity, Amplitude, Modulus, RealAmplitude, RealAmp, ComplexAmplitude, ComplexAmp, LatticeField, LF
export subfield, wrap, square

export natlat, padout, sft, isft
export RealPhase, Generic, Phase, FieldVal, ComplexPhase



#region -----------------------natlat------------------------
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
#endregion

#region -------------------padout, lattice, then array, then field-------------------
"""
    padout(L::AbstractRange, p::Tuple{Int,Int}) 
    Pads an AbstractRange by a certain number of points before the first element and after the last element.
    The padding amount is specified by the tuple `p`, where `p[1]` is the padding before and `p[2]` is the padding after.
    It returns a new range with the specified padding added, scaled by the step of the original range,
    and shifted to start at the correct value. Several overloads are provided for convenience."""
function padout(L::AbstractRange, p::Tuple{Int,Int})
    return (0:sum(p)+length(L)-1) .* step(L) .+ (L[1] - step(L) * p[1])
end
padout(L::AbstractRange, p::Int) = padout(L, (p, p))
padout(L::Lattice{N}, p::NTuple{M,Int}) where {N,M} = M == 2N ? ((padout(L[i], (p[2i-1], p[2i])) for i = 1:N)...,) : error("Bad tuple length.")
padout(L::Lattice{N}, p::NTuple{N,Int}) where {N} = ((padout(L[i], p[i]) for i = 1:N)...,)
padout(L::Lattice{N}, p::Int) where {N} = ((padout(L[i], p) for i = 1:N)...,)

"""
    padout(A::Array{T,N}, p::NTuple{M,Int}, filler=zero(T))
    Pads an N-dimensional array `A` with a certain padding specified by the tuple `p` and uses `filler` as the padding value.
    The `filler` defaults to zero of the same type as the array. Several overloads are provided for convenience.
    """
function padout(A::Array{T,N}, p::NTuple{M,Int}, filler=zero(T)) where {T,N,M}
    M == 2N || error("Bad tuple length.")
    output = fill(filler, size(A) .+ ((p[2i-1] + p[2i] for i = 1:N)...,))
    output[((p[2i-1]+1:size(output, i)-p[2i]) for i = 1:N)...] .= A
    return output
end
padout(A::Array{T,N}, p::NTuple{N,Int}, filler=zero(T)) where {T,N} = padout(A, ((p[ceil(Int, i / 2)] for i = 1:2N)...,), filler)
padout(A::Array{T,N}, p::Int, filler=zero(T)) where {T,N} = padout(A, ((p for i = 1:N)...,), filler)

"""
    padout(f::LF{S,T,N}, p::NTuple{N,Int}, filler=zero(T)) where {S<:FieldVal,T,N}

    Pads a `LatticeField` with specified padding and a filler value. This function extends the `data` field of 
    a `LatticeField` in each dimension according to the padding specified by `p`, and fills the new elements 
    with `filler`, which defaults to the zero value of type `T`.

    # Parameters
    - `f`: An instance of `LatticeField` to be padded.
    - `p`: An N-tuple specifying the padding for each dimension. Each element of the tuple represents the total 
    padding (both before and after) for the corresponding dimension of the lattice.
    - `filler`: An optional value used to fill the new padded elements. Defaults to `zero(T)`.

    # Returns
    A new `LatticeField` instance with the data field padded as specified and the lattice adjusted accordingly.

    This function is particularly useful for operations that require modifying the size of the `LatticeField` data 
    while maintaining its structural integrity, such as in image processing or spatial analysis tasks in optical 
    trapping applications.
    """
padout(f::LF{S,T,N}, p::NTuple{N,Int}, filler=zero(T)) where {S<:FieldVal,T,N} = LF{S,T,N}(padout(f.data, p, filler), padout(f.L, p), f.flambda)

#endregion

#region -------------------sft, isft-------------------
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

#endregion



end # module Core