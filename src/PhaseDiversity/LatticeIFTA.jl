module LatticeIFTA

using FFTW
using ..LatticeTools
using ..Misc
export phasor, gs, gsIter, pdgs, pdgsIter, oneShot

"""
    phasor(z::ComplexF64) -> ComplexF64

    Compute the phasor (unit vector in the complex plane) of a given complex number `z`. The phasor is calculated as `z` divided by its absolute value, which normalizes `z` to have a magnitude of 1. If `z` is zero, the function returns 1.0 (represented as a complex number).

    This function is useful in applications involving complex numbers where only the direction (phase) of the number is relevant, and its magnitude can be normalized.

    Arguments:
    - `z::ComplexF64`: A complex number.

    Returns:
    - `ComplexF64`: The phasor of `z`, which is a complex number with the same phase as `z` but with a magnitude of 1. Returns `1.0 + 0.0im` if `z` is zero.

    """
phasor(z::ComplexF64) = iszero(z) ? one(ComplexF64) : z / abs(z)

#region -------------------Gerchberg-Saxton------------------------------
"""
    gs(U::LF{Modulus,<:Real,N}, V::LF{Modulus,<:Real,N}, nit::Integer, Φ0::LF{ComplexPhase,<:Complex,N}) where N -> LF{ComplexPhase}

    Perform the Gerchberg-Saxton (GS) algorithm iterations on two distributions specified by `U` and `V` for `nit` number of iterations. The `Φ0` parameter is used as the starting phase guess. The function first checks for lattice compatibility between `U`, `V`, and `Φ0`.

    The GS algorithm is used in optics and image processing to iteratively find a phase distribution that, when Fourier transformed, matches a given amplitude distribution.

    Arguments:
    - `U::LF{Modulus,<:Real,N}`: A `LatticeField` representing the first amplitude distribution.
    - `V::LF{Modulus,<:Real,N}`: A `LatticeField` representing the second amplitude distribution.
    - `nit::Integer`: The number of iterations to run the GS algorithm.
    - `Φ0::LF{ComplexPhase,<:Complex,N}`: A `LatticeField` representing the initial phase guess.

    Returns:
    - `LF{ComplexPhase}`: A `LatticeField` representing the complex phase after `nit` iterations of the GS algorithm.

    This function is essential in simulations and calculations where phase retrieval or phase distribution optimization is required.
    """
function gs(U::LF{Modulus,<:Real,N}, V::LF{Modulus,<:Real,N}, nit::Integer, Φ0::LF{ComplexPhase,<:Complex,N}) where {N}
    # Runs nit number of iterations of Gerchberg-Saxton on the distributions specified by U and V.  Φ0 is a starting phase guess.
    ldq(U, V), elq(U, Φ0)                              # Check for compatibility of the lattices
    ft = plan_fft(zeros(ComplexF64, size(U)))
    ift = plan_ifft(zeros(ComplexF64, size(V)))

    guess = ifftshift(U.data .* phasor.(Φ0.data))
    ushift = ifftshift(U.data)
    vshift = ifftshift(V.data)
    for i = 1:nit
        guess = gsIter(guess, ushift, vshift, ft, ift)
    end
    return LF{ComplexPhase}(fftshift(phasor.(guess)), U.L)
end

"""
    gsIter(guess::Array{Complex{Float64},N}, u::Array{Float64, N}, v::Array{Float64,N}, ft::FFTW.cFFTWPlan, ift::AbstractFFTs.ScaledPlan) where N

    Perform a single iteration of the Gerchberg-Saxton algorithm. This function applies Fourier transforms to the `guess` array, modifies it based on the amplitude constraints `u` and `v`, and then applies the inverse Fourier transform.

    Arguments:
    - `guess`: The current guess for the complex field.
    - `u`: The target amplitude in the object plane.
    - `v`: The target amplitude in the Fourier plane.
    - `ft`: The FFT plan for forward transformation.
    - `ift`: The FFT plan for inverse transformation.
    """
function gsIter(guess::Array{Complex{Float64},N}, u::Array{Float64,N}, v::Array{Float64,N}, ft::FFTW.cFFTWPlan, ift::AbstractFFTs.ScaledPlan) where {N}
    # Single iteration of Gerchberg-Saxton using provided Fourier operators
    return u .* phasor.(ift * (phasor.(ft * guess) .* v))
end

"""
    gs(U::LF{Modulus,<:Real,N}, V::LF{Modulus,<:Real,N}, nit::Integer, Φ0::LF{RealPhase,<:Real,N}) where N

    Gerchberg-Saxton algorithm for `LatticeField` with `RealPhase`. It wraps the `RealPhase` to `ComplexPhase` before running the algorithm.

    Arguments:
    - `U`, `V`: `LatticeField` instances representing amplitude distributions.
    - `nit`: Number of iterations.
    - `Φ0`: Initial phase guess as `RealPhase`.
    """
function gs(U::LF{Modulus,<:Real,N}, V::LF{Modulus,<:Real,N}, nit::Integer, Φ0::LF{RealPhase,<:Real,N}) where {N}
    return gs(U, V, nit, wrap(Φ0))
end

"""
    gs(U::LF{Modulus,<:Real,N}, V::LF{Modulus,<:Real,N}, nit::Integer) where N

    Gerchberg-Saxton algorithm for `LatticeField` without an initial phase guess. It initializes a random phase guess.

    Arguments:
    - `U`, `V`: `LatticeField` instances representing amplitude distributions.
    - `nit`: Number of iterations.
    """
function gs(U::LF{Modulus,<:Real,N}, V::LF{Modulus,<:Real,N}, nit::Integer) where {N}
    return gs(U, V, nit, LF{RealPhase}(rand(size(U)...), U.L, U.flambda))
end

"""
    gs(U::LF{Intensity,<:Real,N}, V::LF{Intensity,<:Real,N}, nit::Integer, Φ0::Union{Nothing,LF{ComplexPhase,<:Complex,N},LF{RealPhase,<:Real,N}}=nothing) where N

    Gerchberg-Saxton algorithm for `LatticeField` with `Intensity`. It converts intensity to modulus by taking the square root before running the algorithm.

    Arguments:
    - `U`, `V`: `LatticeField` instances representing intensity distributions.
    - `nit`: Number of iterations.
    - `Φ0`: Optional initial phase guess.
    """
function gs(U::LF{Intensity,<:Real,N}, V::LF{Intensity,<:Real,N}, nit::Integer, Φ0::Union{Nothing,LF{ComplexPhase,<:Complex,N},LF{RealPhase,<:Real,N}}=nothing) where {N}
    return gs(sqrt(U), sqrt(V), nit, Φ0)
end
#endregion

#region ----------------Phase Diversity --------------------


"""
    pdgs(imgs::NTuple{M,LF{Modulus,<:Number,N}}, divPhases::NTuple{M,LF{<:Phase,<:Number,N}}, nit::Integer, beamGuess::LF{ComplexAmp,<:Complex,N}) where {M,N}

    Perform the phase diversity Gerchberg-Saxton (PDGS) algorithm. This function is designed for a set of images (`imgs`) with corresponding diversity phases (`divPhases`). It iteratively refines an initial guess of the beam (`beamGuess`) over `nit` iterations.

    Arguments:
    - `imgs`: A tuple of `LatticeField` instances representing the modulus of the beam distributions.
    - `divPhases`: A tuple of `LatticeField` instances representing phase diversities.
    - `nit`: Number of iterations.
    - `beamGuess`: Initial guess for the complex amplitude of the beam.

    Returns:
    - `LF{ComplexAmp}`: Refined complex amplitude of the beam after PDGS iterations.
    """
function pdgs(imgs::NTuple{M,LF{Modulus,<:Number,N}}, divPhases::NTuple{M,LF{<:Phase,<:Number,N}}, nit::Integer, beamGuess::LF{ComplexAmp,<:Complex,N}) where {M,N}
    s = size(imgs[1])
    all(size(i) == s for i in imgs) && all(size(i) == size(beamGuess) for i in divPhases) || error("Input size mismatch.")
    [elq(d, beamGuess) for d in divPhases], [ldq(i, beamGuess) for i in imgs]
    ft = plan_fft(zeros(s))
    ift = plan_ifft(zeros(s))
    guess = ifftshift(beamGuess.data)
    phis = ((ifftshift(wrap(divPhases[i]).data .* dualPhase(imgs[i].L, imgs[i].flambda).data) * 1.0 for i = 1:M)...,)::NTuple{M,Array{ComplexF64,N}}
    mods = ((ifftshift(i.data) * 1.0 for i in imgs)...,)::NTuple{M,Array{Float64,N}}
    for j = 1:nit
        guess = pdgsIter(guess, phis, mods, ft, ift)
    end
    return LF{ComplexAmp}(fftshift(guess), beamGuess.L, beamGuess.flambda)
end

"""
    pdgs(imgs::NTuple{M,LF{Intensity,<:Number,N}}, divPhases::NTuple{M,LF{<:Phase,<:Number,N}}, nit::Integer, beamGuess::LF{ComplexAmp,<:Complex,N}) where {M,N}

    Variation of the PDGS algorithm for `LatticeField` instances with `Intensity`. It converts the intensity to modulus by taking the square root before running the PDGS algorithm.

    Arguments:
    - `imgs`: A tuple of `LatticeField` instances representing the intensity of the beam distributions.
    - `divPhases`: A tuple of `LatticeField` instances representing phase diversities.
    - `nit`: Number of iterations.
    - `beamGuess`: Initial guess for the complex amplitude of the beam.

    Returns:
    - `LF{ComplexAmp}`: Refined complex amplitude of the beam after PDGS iterations.
    """
function pdgs(imgs::NTuple{M,LF{Intensity,<:Number,N}}, divPhases::NTuple{M,LF{<:Phase,<:Number,N}}, nit::Integer, beamGuess::LF{ComplexAmp,<:Complex,N}) where {M,N}
    return pdgs(Tuple(sqrt(i) for i in imgs), divPhases, nit, beamGuess)
end

"""
    pdgsIter(guess::Array{Complex{Float64},N}, phis::NTuple{M,Array{ComplexF64,N}}, mods::NTuple{M,Array{Float64,N}}, ft::FFTW.cFFTWPlan, ift::AbstractFFTs.ScaledPlan) where {M,N}

    Perform a single iteration of the phase diversity Gerchberg-Saxton algorithm. This function updates the `guess` of the beam's complex amplitude based on the modulus constraints (`mods`) and phase diversities (`phis`).

    Arguments:
    - `guess`: The current guess of the beam's complex amplitude.
    - `phis`: Tuple of phase diversity arrays.
    - `mods`: Tuple of modulus arrays.
    - `ft`: FFT plan for forward transformation.
    - `ift`: FFT plan for inverse transformation.
    """
function pdgsIter(guess::Array{Complex{Float64},N}, phis::NTuple{M,Array{ComplexF64,N}}, mods::NTuple{M,Array{Float64,N}},
    ft::FFTW.cFFTWPlan, ift::AbstractFFTs.ScaledPlan) where {M,N}
    # Single iteration of Gerchberg-Saxton phase diversity using provided Fourier operators
    return +((ift * (mods[i] .* phasor.(ft * (guess .* phis[i]))) .* conj.(phis[i]) for i = 1:M)...) ./ M
end

#endregion


"""
    oneShot(img::LF{Intensity,<:Real,N}, alpha::Real, beta::NTuple{N,Real}, sc::Real) where N

    Perform a one-shot phase retrieval process on an image `img` using specified parameters `alpha`, `beta`, and `sc`. This function is designed for applications in optics and image processing where quick and efficient phase retrieval from a single intensity image is required.

    Arguments:
    - `img`: A `LatticeField` representing the intensity distribution of the image.
    - `alpha`: A real number parameter used in the phase retrieval process.
    - `beta`: An `N`-tuple of real numbers representing additional parameters for each dimension.
    - `sc`: A scaling factor used in the phase retrieval process.

    Returns:
    - `LF{ComplexAmp}`: A `LatticeField` representing the complex amplitude after the phase retrieval process.

    The function works by first calculating a dual lattice and then applying a series of transformations to retrieve the phase information from the intensity image. The result is a complex amplitude field that approximates the original field of the image.
    """
function oneShot(img::LF{Intensity,<:Real,N}, alpha::Real, beta::NTuple{N,Real}, sc::Real) where {N}
    dL = dualShiftLattice(img.L, img.flambda)
    xc = beta .* (img.flambda / sc)
    divPhase = LF{RealPhase}([alpha / (2 * sc^2) * sum(dL[j][I[j]]^2 for j = 1:N) + sum(beta[j] * dL[j][I[j]] for j = 1:N) / sc for I in CartesianIndices(size(b))])
    dualDivPhase = LF{RealPhase}([sc^2 / (2 * alpha * flambda^2) * sum((L[j][I[j]] - xc[j])^2 for j = 1:N) for I in CartesianIndices(size(img))], img.L, img.flambda)
    x = sqrt(img) * dualDivPhase
    return LF{ComplexAmp}(x.data |> ifftshift |> ifft |> fftshift, dL, x.flambda) * conj(divPhase)
end


end # module LatticeIFTA