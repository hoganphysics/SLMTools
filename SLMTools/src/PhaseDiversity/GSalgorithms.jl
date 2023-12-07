# Gerchberg-Saxton algorithms for phase retrieval.
using FileIO, Images, Plots, Interpolations, FFTW, Statistics
"""
    phasor(z::Complex)

Converts a complex number `z` to its phasor representation. If `z` is zero, it returns zero; otherwise, it returns `z` normalized by its absolute value.

# Arguments
- `z::Complex`: A complex number.

# Returns
- `Complex`: The phasor representation of `z`, of unit length.
"""
phasor(z::Complex) = iszero(z) ? zero(z * 1.0) : z * 1.0 / abs(z)

"""
    arraySizeChecks(arrs...)

    Performs size consistency checks on a variable number of arguments `arrs`. Each argument must be either a 2D array (matrix) or a vector of 2D arrays. The function checks that all matrices have the same size and, if vectors are provided, that all vectors have the same length.

    # Arguments
    - `arrs...`: A variable number of arrays or vectors of arrays.

    # Returns
    - `(Tuple{Union{Nothing, Int}, Union{Nothing, Int}})`: A tuple containing the size of the matrices and the length of the lists, if applicable. If one is not defined, it returns `nothing` in its place.

    # Errors
    - Throws an error if the input arguments do not meet the specified criteria.
"""
function arraySizeChecks(arrs...)
    # Checks that all inputs are either matrices or lists of matrices, where each list has the same length and all matrices have the same size. 
    # Returns (size of matrices, length of lists) if both are defined.  If one is not defined, it returns `nothing` in its place. 
    all([(a isa Array{R,2} where {R<:Number}) || (a isa Vector{Array{R,2}} where {R<:Number}) for a in arrs]) || error("All arguments must be matrices or vectors of matrices.")
    vs = [a for a in arrs if (a isa Vector)]
    ms = [a for a in arrs if (a isa Matrix)]
    (length(vs) + length(ms) == length(arrs)) || error("Something went wrong")
    n, s = nothing, nothing
    if length(vs) > 0
        n = length(vs[1])
        s = size(vs[1][1])
        all(length(v) == n for v in vs) || error("All lists must have equal length.")
    end
    if length(ms) > 0
        s = size(ms[1])
    end
    all(size(v[1]) == s for v in vs) && all(size(m) == s for m in ms) || error("All matrices must have equal size.")
    return s, n
end

"""
    ampSwap(ampArray, phaseArray)

Performs element-wise multiplication of the amplitudes from `ampArray` with the phasors computed from `phaseArray`.

# Arguments
- `ampArray`: An array containing amplitude values.
- `phaseArray`: An array containing phase values.

# Returns
- `Array`: An array of the product of amplitudes and corresponding phasors.
"""
function ampSwap(ampArray, phaseArray)
    return abs.(ampArray) .* phasor.(phaseArray)
end

"""
    pdgs(imgs, phases, nit, guess0)

Performs a Gerchberg-Saxton style phase diversity inversion algorithm. It iteratively updates an initial guess `guess0` using a set of images `imgs` and corresponding phase masks `phases`.

# Arguments
- `imgs`: An array of images.
- `phases`: An array of corresponding phase masks.
- `nit`: Number of iterations to perform.
- `guess0`: Initial guess for the iterative process.

# Returns
- `(Array, Array)`: The final guess after iterations and an array of variance values for each iteration.

# Errors
- Throws an error if input arrays do not pass size checks or if phase masks are not complex.
"""
function pdgs(imgs, phases, nit, guess0)
    # Error checking
    sz, ni = arraySizeChecks(imgs, phases, guess0)            # sz, ni = (size of image, number of images)
    isa(phases[1][1, 1], Complex) || error("Phases should be complex.")

    ft = plan_fft(zeros(size(imgs[1])))
    ift = plan_ifft(zeros(size(imgs[1])))
    guess = guess0
    vars = zeros(nit)
    amps = [sqrt.(img / sum(img)) for img in imgs]
    for i = 1:nit
        updates = [(ift * ifftshift(ampSwap(amps[j], fftshift(ft * (guess .* phases[j]))))) .* conj.(phases[j]) for j = 1:ni]
        updates ./= [sqrt.(sum(abs.(u) .^ 2)) for u in updates]    # Normalize to unit total power for each image.
        guess = mean(updates)
        vars[i] = ([abs.(u - guess) .^ 2 for u in updates] |> sum |> sum) / (ni * prod(sz))
    end
    return guess, vars
end

"""
    gs(ampin, ampout, nit, guess0; divPhase=nothing)

Implements the standard Gerchberg-Saxton algorithm with optional diversity phase compensation. It iteratively refines a guess `guess0` to match the amplitude constraints of `ampin` and `ampout`.

# Arguments
- `ampin`: Input amplitude constraints.
- `ampout`: Output amplitude constraints.
- `nit`: Number of iterations to perform.
- `guess0`: Initial guess for the iterative process.
- `divPhase`: Optional diversity phase for compensation.

# Returns
- `(Array, Array)`: The final refined guess and an array of variance values for each iteration.

# Notes
- If `divPhase` is provided, it is used for diversity phase compensation.
"""
function gs(ampin, ampout, nit, guess0; divPhase=nothing)
    ft = plan_fft(zeros(size(ampin)))
    ift = plan_ifft(zeros(size(ampout)))
    if isnothing(divPhase)
        guess = ifftshift(guess0) ./ sqrt(sum(abs.(guess0) .^ 2))
    else
        guess = ifftshift(guess0 .* divPhase) ./ sqrt(sum(abs.(guess0) .^ 2))
    end
    aos = ifftshift(abs.(ampout)) ./ sqrt(sum(abs.(ampout) .^ 2))
    ais = ifftshift(ampin) ./ sqrt(sum(abs.(ampin) .^ 2))
    vars = zeros(nit)
    for i = 1:nit
        ftguess = ft * guess
        ftguess ./= sqrt(sum(abs.(ftguess) .^ 2))
        guess = ampSwap(ais, ift * ampSwap(aos, ftguess))
        vars[i] = sum(abs.(aos - abs.(ftguess)) .^ 2)
    end
    if !isnothing(divPhase)
        guess .*= conj.(ifftshift(divPhase))
    end
    return fftshift(guess), vars
end

"""
    croppedPDGS(imgs::Vector{Array{T,N}}, lattices::Vector{S}, phases::Vector{Array{R,N}}, nit::Integer, guess0::Union{Array{R,N},Nothing}) where {T<:Real, S<:Lattice{N}, R<:Complex, N}

Implements a cropped Phase Diversity Gerchberg-Saxton (GS) algorithm. This function takes a set of images, corresponding lattices, and phases, and iteratively refines an initial guess to estimate the wavefront or phase distribution. It is particularly useful in adaptive optics and imaging through turbulent media.

# Arguments
- `imgs`: A vector of images (arrays) used in the phase diversity estimation.
- `lattices`: A vector of lattice structures corresponding to each image.
- `phases`: A vector of initial phase estimates corresponding to each image.
- `nit`: The number of iterations to perform in the GS algorithm.
- `guess0`: An initial guess for the wavefront estimate. If `Nothing`, a zero array is used.

# Returns
- An array representing the estimated wavefront after `nit` iterations of the GS algorithm.

# Notes
- This function assumes that the lengths and sizes of `imgs`, `lattices`, and `phases` are consistent. It throws an error if there is a mismatch in lengths or sizes.
- The algorithm uses Fast Fourier Transforms (FFTs) for its computations.
"""
function croppedPDGS(imgs::Vector{Array{T,N}}, lattices::Vector{S}, phases::Vector{Array{R,N}}, nit::Integer, guess0::Union{Array{R,N},Nothing}) where {T<:Real,S<:Lattice{N},R<:Complex} where {N}
    m = length(imgs)
    s = size(imgs[1])
    length(imgs) == length(lattices) == length(phases) || error("Input length mismatch.")
    all(size(i) == s for i in imgs) && all(length.(l) == s for l in lattices) && all(size(p) == s for p in phases) || error("Input size mismatch.")
    isnothing(guess0) && (guess0 = zeros(s))

    ft = plan_fft(zeros(s))
    ift = plan_ifft(zeros(s))
    guess = ifftshift(guess0)
    ϕs = ifftshift.([phases[j] .* dualPhase(lattices[j]) for j = 1:m])
    amps = ifftshift.([sqrt.(i ./ sum(i)) for i in imgs])

    for i = 1:nit
        θs = [exp.(im * angle.(ft * (guess .* ϕs[j]))) for j = 1:m]
        guess = sum([(ift * (θs[j] .* amps[j])) .* conj.(ϕs[j]) for j = 1:m]) ./ m
    end
    return fftshift(guess)
end

"""
    pdError1(beam::Array{ComplexF64,N}, pdis::Vector{Array{T,N}}, ψs::Vector{Array{ComplexF64,N}}) where {N, T<:Real}

Calculates the phase diversity error for a given beam and set of phase diversity images (`pdis`). This error metric is used to evaluate the quality of phase recovery in phase diversity techniques.

# Arguments
- `beam`: A complex array representing the beam.
- `pdis`: A vector of arrays representing phase diversity images.
- `ψs`: A vector of complex arrays representing phase adjustments.

# Returns
- A scalar value representing the phase diversity error.
"""
function pdError1(beam::Array{ComplexF64,N}, pdis::Vector{Array{T,N}}, ψs::Vector{Array{ComplexF64,N}}) where {N,T<:Real}
    sqrt(sum(sum((nabs(sft(beam .* ψs[j])) .^ 2 - (pdis[j] ./ sum(pdis[j]))) .^ 2) for j = 1:length(pdis)))
end

"""
    cpdError1(beam::Array{ComplexF64,N}, pdis::Vector{Array{T,N}}, ψs::Vector{Array{ComplexF64,N}}, lattices::Vector{L}) where {N, T<:Real, L<:Lattice}

Calculates a corrected phase diversity error for a given beam, phase diversity images, phase adjustments, and lattice structures. This function extends the `pdError1` by incorporating the effect of the lattice structures.

# Arguments
- `beam`: A complex array representing the beam.
- `pdis`: A vector of arrays representing phase diversity images.
- `ψs`: A vector of complex arrays representing phase adjustments.
- `lattices`: A vector of lattice structures.

# Returns
- A scalar value representing the corrected phase diversity error.
"""
function cpdError1(beam::Array{ComplexF64,N}, pdis::Vector{Array{T,N}}, ψs::Vector{Array{ComplexF64,N}}, lattices::Vector{L}) where {N,T<:Real,L<:Lattice}
    sqrt(sum(sum((nabs(sft(beam .* ψs[j] .* dualPhase(lattices[j]))) .^ 2 - (pdis[j] ./ sum(pdis[j]))) .^ 2) for j in 1:length(pdis)))
end

"""
    gsRefine(μ::AbstractArray{<:Number,N}, ν::AbstractArray{<:Number,N}, nit::Integer, Φ::AbstractArray{ComplexF64,N}) where N

Refine an initial guess of a phase function using the Gerchberg-Saxton algorithm.

# Arguments
- `μ::AbstractArray{<:Number,N}`: A source distribution array.
- `ν::AbstractArray{<:Number,N}`: A target distribution array.
- `nit::Integer`: The number of iterations for the algorithm.
- `Φ::AbstractArray{ComplexF64,N}`: The initial guess for the phase function.

# Returns
- An array representing the refined phase function.

This function takes an initial guess for the phase function and iteratively refines it using the Gerchberg-Saxton algorithm to better approximate the phase that transforms the source distribution `μ` into the target distribution `ν`.
"""
function gsRefine(μ::AbstractArray{<:Number,N}, ν::AbstractArray{<:Number,N}, nit::Integer, Φ::AbstractArray{ComplexF64,N}) where {N}
    guess0 = sqrt.(μ) .* Φ
    ΦR, vs = gs(sqrt.(μ), sqrt.(ν), nit, guess0)
    return ΦR
end