#=
Utility functions for working with LatticeFields.
Requires:
    using FFTW: fftshift, ifftshift, fft
    using LinearAlgebra 
    using ToeplitzMatrices 

Notes: 
    All the dependencies are just for Andrii's Wigner and FrFT functions. 
=#

export Lattice, elq, RealPhase, Generic, Phase, FieldVal, ComplexPhase, UPhase, UnwrappedPhase, S1Phase
export Intensity, Amplitude, Modulus, RealAmplitude, RealAmp, ComplexAmplitude, ComplexAmp, LatticeField, LF
export subfield, wrap, square, sublattice, normalizeLF, phasor
export natrange, natlat, padout, latticeDisplacement, toDim, naturalize, r2, ldot, Nyquist




#-----------------------natlat------------------------
""" 
    natlat(n::Int) 
    Create a "natural" 1D lattice for a given integer n, i.e. a lattice which is self-dual under DFT.
	"""
function natrange(n::Int)
    return (-floor(n / 2):floor((n - 1) / 2)) ./ sqrt(n)
end

"""
    natlat(ns::NTuple{N,Integer}) 
    Make a natlat for an N-dimensional lattice specified by a tuple of integers."""
function natlat(ns::NTuple{N,Integer}) where {N}
    return ((natrange(n) for n in ns)...,)
end

"""
    natlat(ns::Integer...) 
    Make a natlat for an N-dimensional lattice specified by N integer arguments."""
function natlat(ns::Real...)
    return natlat(ns)
end

"""
    naturalize(f::LF{S,T,N}) where {S<:FieldVal, T, N}

    Convert a LatticeField to be defined on its natural lattice.

    This function takes a `LF{S,T,N}` object and converts it to be defined on the natural lattice that corresponds to its size.

    # Arguments
    - `f::LF{S,T,N}`: A LatticeField object to be naturalized.

    # Returns
    - A `LF{S}` object with the same data as `f`, but defined on its natural lattice.
    """
function naturalize(f::LF{S,T,N}) where {S<:FieldVal,T,N}
    # Converts LatticeField to be difined on the natural lattice of the same size. 
    return LF{S}(f.data, natlat(size(f)), 1.0)
end


#-------------------padout, lattice, then array, then field-------------------
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
    """
padout(f::LF{S,T,N}, p::NTuple{N,Int}, filler=zero(T)) where {S<:FieldVal,T,N} = LF{S,T,N}(padout(f.data, p, filler), padout(f.L, p), f.flambda)

#-------------------lattice displacement and toDim-------------------

"""
    latticeDisplacement(L::Lattice)

    Calculate the displacement for each dimension of a lattice `L`. This function determines the offset from the origin for each dimension of the lattice.

    Arguments:
    - `L::Lattice`: The lattice to calculate the displacement for.

    Returns:
    - Array: Displacements for each dimension of the lattice.

    Usage:
    displacements = latticeDisplacement(myLattice)
    """
function latticeDisplacement(L::Lattice)
    return [l[1] + floor(length(l) / 2) * step(l) for l in L]
end


"""
    toDim(v, d::Int, n::Int)

    Reshape a vector `v` into an `n`-dimensional array, with non-singleton dimension at `d`. This is useful for applying operations intended for one dimension across a specific dimension in a higher-dimensional space.

    Arguments:
    - `v`: The vector to be reshaped.
    - `d::Int`: The dimension which `v` should span.
    - `n::Int`: The number of dimensions of the target array.

    Returns:
    - Array: The reshaped array.

    Usage:
    reshapedArray = toDim(vector, dimension, totalDimensions)
    """
function toDim(v, d::Int, n::Int)
    reshape(v, ((i == d ? length(v) : 1) for i = 1:n)...)
end

#------------------- Lattice vector math -------------------

"""
	r2(x::Lattice{N}) 
	
	returns an array of the L2 norm squared for each point associated with a lattice. 
	"""
r2(x::Lattice{N}) where N = .+( (toDim(x[i].^2,i,N) for i=1:N)... )

"""
	ldot(v,x::Lattice{N})
	
	Computes the dot product of a vector with the coordinates associated with a lattice.  
	Has methods for all likely combinations of vectors and lattices, i.e. either order and 
	with the vector being represented by a Vector or a Tuple. 
	"""
ldot(v::NTuple{N,Real},x::Lattice{N}) where N = .+( (toDim(x[i] * v[i],i,N) for i=1:N)... )
ldot(x::Lattice{N},v::NTuple{N,Real}) where N = ldot(v,x)
ldot(v::Vector{T},x::Lattice{N}) where {T,N} = (length(v)==N || error("Vector length != Lattice dimension."); .+( (toDim(x[i] * v[i],i,N) for i=1:N)... ))
ldot(x::Lattice{N},v::Vector{T}) where {N,T} = ldot(v,x)

#------------------- Nyquist frequency -------------------

"""
	Nyquist(L::Lattice{N}) 
	
	Returns a tuple of the Nyquist frequencies of the lattice L in each dimension.  Useful for making phase ramps for the SLM; 
	often you'll want the phase ramp to be Nyquist(L) .* (1/2,0).
	"""
function Nyquist(L::Lattice{N}) where N
    return 1 ./ (2 .* step.(L))
end
#------------------- Basis transformations -------------------

"""
    shiftedDFTBasis(n::Integer)

    Generate a shifted Discrete Fourier Transform (DFT) basis matrix of size n×n.
    The basis is centered around the origin, with the zero frequency at the center.

    # Arguments
    - `n::Integer`: Size of the DFT matrix

    # Returns
    - `W::Matrix{Complex}`: The shifted DFT basis matrix
"""
function shiftedDFTBasis(n::Integer)
    shift = n ÷ 2
    W = [1/sqrt(n)*exp(-1im * 2 * π * (j-shift)*(k-shift)/n) for j in 0:n-1, k in 0:n-1];
    return W
end 

"""
    hermiteBasis(n::Integer)

    Generate a discrete Hermite-Gaussian basis, which is simultanously an eigenbasis
    of the shiftedDFTBasis.

    # Arguments
    - `n::Integer`: Size of the basis set

    # Returns
    - `H::Matrix`: The Hermite-Gaussian basis matrix
    - `λ::Vector`: The eigenvalues of the quantum harmonic oscillator
"""
function hermiteBasis(n::Integer)
    # scale L by 2pi (eigenvalue convention)
    L = natlat(n)[1]
    L = L .* 2π

    # Generate a shifted DFT matrix
    W = shiftedDFTBasis(L)
    
    # Form a Hamiltonian matrix
    Q = diagm(L)
    P = W*Q*W'
    QHO = P'*P + Q'*Q

    # Perform the diagonalization
    λ, H1 = eigen(QHO);

    # QR decomposition to make them real-valued
    H1 = real((H1 + conj(H1))/2)
    H,_ = qr(H1)
    H = Matrix(H)
    
    return H, λ
end

"""
    FrFTBasis(n::Integer, α::Real)

    Generate a Fractional Fourier Transform (FrFT) basis matrix.
    The FrFT is a generalization of the Fourier transform with fractional order α.

    # Arguments
    - `n::Integer`: Size of the FrFT matrix
    - `α::Real`: Fractional order of the transform (in units of π/2)

    # Returns
    - `Wα::Matrix`: The FrFT basis matrix
"""
function FrFTBasis(n::Integer, α::Real)
    # generates hermite basis
    W = shiftedDFTBasis(n)
    H, λ = hermiteBasis(n)
    
    # generates FrFT matrix
    λw = diag(H' * W * H)
    Λ = diagm(λw.^(2α/π))
    Wα = H * Λ * transpose(H)
    return Wα
end

"""
    wigner_fft(f::Vector{T}) where T<:Number

    Compute the Wigner distribution function using FFT method.
    The Wigner distribution has the first coordinate that lives on the 
    same lattice as the input data. The second coordinate is on the lattice 
    that has twice as many points, and is half as big in the extent.

    # Arguments
    - `f::Vector{T}`: Input signal vector

    # Returns
    - `w::Matrix`: The Wigner distribution matrix
"""
function wigner_fft(f::Vector{T}) where T<:Number
    # Reshape function
    f = reshape(f,1,:).|>ComplexF64
    n = length(f[:])
    
    # Construct r1 and r2
    r1 = hcat([0], reverse(conj.(f)), zeros(1, n - 1))
    r2 = hcat([0], f, zeros(1, n - 1))

    # Create Toeplitz matrices
    w = Toeplitz(zeros(n), r1[:]) .* reverse(Toeplitz(zeros(n), r2[:]), dims=1)
    
    # compute FT 
    w = 1/(2n) * ifftshift(fft(fftshift(w, 2),2),2)
    
    return real(w)
end
