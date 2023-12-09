"""
    getMoment(A::AbstractArray{<:Number,N}, m::NTuple{M,Integer}) where {N,M}

    Compute the moment of a given array `A`.

    # Arguments
    - `A::AbstractArray{<:Number,N}`: A multidimensional array representing a discrete function or data set.
    - `m::NTuple{M,Integer}`: A tuple specifying the indices for the moment calculation. For instance, `(i_1, i_2, ..., i_M)` indicates that the moment `<x_{i_1} x_{i_2} ... x_{i_M}>` should be calculated.

    # Returns
    - The calculated moment of the array `A` based on the specified indices in `m`.

    # Details
    This function calculates the moment of the array `A` along the specified dimensions in `m`. It internally handles non-dimensional and dimensional aspects of the array to focus the moment calculation on the relevant axes.
    """
function getMoment(A::AbstractArray{<:Number,N}, m::NTuple{M,Integer}) where {N,M}
    # Returns the moment <x_{i_1} x_{i_2} ... x_{i_M}> of A, where m = (i_1,i_2,...,i_M)
    pows = zeros(Int, N)
    for j in m
        pows[j] += 1
    end
    nondims = filter(x -> pows[x] == 0, 1:N)
    indims = filter(x -> pows[x] != 0, 1:N)
    B = sum(A, dims=nondims)
    return sum(B .* [prod([I[j]^pows[j] for j in indims]) for I in CartesianIndices(size(B))])
end


#TODO do i need this?
"""
    downsample(μ::AbstractArray{T,N}, r::Int; reducer=sum) where {T, N}

    Downsampling operation on an N-dimensional array `μ` by a factor of `r` using a specified reduction function (default is sum). This function reduces the size of the array by grouping and aggregating elements.

    # Arguments
    - `μ`: An N-dimensional array to be downsampled.
    - `r`: Downsampling factor, which must evenly divide each dimension of `μ`.
    - `reducer`: A function to aggregate elements in each group (default is `sum`).

    # Returns
    - A downsampled version of the input array `μ`.

    # Errors
    - Throws an error if `r` does not evenly divide each dimension of `μ`.
"""
function downsample(μ::AbstractArray{T,N}, r::Int; reducer=sum) where {T,N}
    size(μ) .% r == 0 .* size(μ) || error("Downsample factor must divide array size.")
    CR = CartesianIndices(((0:r-1 for i = 1:N)...,))
    CD = CartesianIndices(((1:r:s for s in size(μ))...,))
    return [reducer(μ[I.+CR]) for I in CD]
end

#TODO do i need this?
"""
    downsample(L::Lattice{N}, r::Int) where {N}

    Downsampling operation for a lattice `L` by a factor of `r`. This function reduces the number of points in the lattice by selecting points at regular intervals based on the downsampling factor.

    # Arguments
    - `L`: A lattice to be downsampled.
    - `r`: Downsampling factor, which must evenly divide the length of each dimension of `L`.

    # Returns
    - A downsampled version of the lattice `L`.

    # Errors
    - Throws an error if `r` does not evenly divide the length of each dimension of `L`.
"""
function downsample(L::Lattice{N}, r::Int) where {N}
    length.(L) .% r == 0 .* length.(L) || error("Downsample factor must divide lattice size.")
    return ((((0:r:length(l)-1) .+ (r - 1) / 2) .* step(l) .+ l[1] for l in L)...,)
end




#region ---------------GSAlgorthims.jl code before refactor----------------
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

#endregion



#region -------estimateInverse ---------


"""
    estimateInverseImage(interp, alpha::Real, slmGrid::Lattice, flambda::Real)

Estimate the inverse image from single diversity images using the Fourier transform method.

# Arguments
- `interp::Array`: The intensity of the input beam represented as an array.
- `alpha::Real`: Scaling factor for the quadratic phase term in the diversity phase.
- `slmGrid::Lattice`: Spatial Light Modulator (SLM) grid, representing the x and y coordinates in the lattice.
- `flambda::Real`: Wavelength scaling factor.

# Returns
`Array`: Returns the estimated intensity distribution of the inverse image.

# Explanation
The function applies a quadratic phase term exp(2πiα(x^2 + y^2) / 2) to the square root of the interpolated intensity.
Then, it performs an inverse Fourier transform to estimate the inverse image. The `fftshift` and `ifftshift` functions
are used to center the Fourier and inverse Fourier transforms.
"""
function estimateInverseImage(interp, alpha::Real, slmGrid::Lattice, flambda::Real)
    # We take the diversity phase to be exp(2*pi*im*alpha*(x^2+y^2) / 2)
    dL = dualShiftLattice(slmGrid, flambda)
    b = sqrt.(interp) .* [exp(-2 * pi * im * (x^2 + y^2) / (2 * alpha * flambda^2)) for x in dL[1], y in dL[2]]
    return fftshift(abs.(ifft(ifftshift(b))) .^ 2)
end

"""
    estimateInverseField(interp, alpha::Real, slmGrid::Lattice, flambda::Real)

Estimate the inverse field from single diversity images using Fourier transform method.

# Arguments
- `interp::Array`: The intensity of the input beam represented as an array.
- `alpha::Real`: Scaling factor for the quadratic phase term in the diversity phase.
- `slmGrid::Lattice`: Spatial Light Modulator (SLM) grid, representing the x and y coordinates in the lattice.
- `flambda::Real`: Wavelength scaling factor.

# Returns
`Array`: Returns the estimated field distribution of the inverse field.

# Explanation
This function is similar to `estimateInverseImage` but returns the complex field instead of the intensity.
It applies the diversity phase to the square root of the interpolated intensity, performs an inverse Fourier
transform, and multiplies it with a compensatory exponential phase term.
"""
function estimateInverseField(interp, alpha::Real, slmGrid::Lattice, flambda::Real)
    dL = dualShiftLattice(slmGrid, flambda)
    b = sqrt.(interp) .* [exp(-2 * pi * im * (x^2 + y^2) / (2 * alpha * flambda^2)) for x in dL[1], y in dL[2]]
    return fftshift(ifft(ifftshift(b))) .* ([exp(-2 * pi * im * (alpha / 2) * (x^2 + y^2)) for x in slmGrid[1], y in slmGrid[2]])
end

"""
    estimateInverseField(interp, alpha::Real, beta::Vector, slmGrid::Lattice, flambda::Real, sc::Real)

Advanced estimation of the inverse field considering lateral shift.

# Arguments
- `interp::Array`: The intensity of the input beam represented as an array.
- `alpha::Real`: Scaling factor for the quadratic phase term in the diversity phase.
- `beta::Vector`: Lateral shift parameters as a vector [x_shift, y_shift].
- `slmGrid::Lattice`: Spatial Light Modulator (SLM) grid, representing the x and y coordinates in the lattice.
- `flambda::Real`: Wavelength scaling factor.
- `sc::Real`: Scale factor for beta adjustments.

# Returns
`Array`: Returns the estimated field distribution of the inverse field with lateral shift considerations.

# Explanation
This function extends `estimateInverseField` by adding a lateral shift to the quadratic phase term.
The lateral shift is determined by `beta`, scaled by `flambda` and `sc`. It is useful for cases
where lateral shifts in the beam need to be accounted for in the inverse field estimation.
"""
function estimateInverseField(interp, alpha::Real, beta::Vector, slmGrid::Lattice, flambda::Real, sc::Real)
    dL = dualShiftLattice(slmGrid, flambda)
    x0, y0 = beta * flambda / sc
    b = sqrt.(interp) .* [exp(-2 * pi * im * ((x - x0)^2 + (y - y0)^2) * sc^2 / (2 * alpha * flambda^2)) for x in dL[1], y in dL[2]]
    return fftshift(ifft(ifftshift(b))) .* ([exp(-2 * pi * im * (alpha * (x^2 + y^2) / (2 * sc^2) + (beta[1] * x + beta[2] * y) / sc)) for x in slmGrid[1], y in slmGrid[2]])
end


#endregion