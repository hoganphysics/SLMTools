module LatticeCore
using FFTW: fftshift, ifftshift
export Lattice, natlat, padout, sft, isft

"""
    Lattice{N}

    A type alias for representing an N-dimensional lattice.

    Each element in the `NTuple` is an `AbstractRange`, representing a range in each dimension of the lattice.

    # Type Parameters
    - `N`: The number of dimensions in the lattice.

    # Examples
    # (Omitted as per your request)
    """
const Lattice{N} = NTuple{N,AbstractRange} where {N}


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