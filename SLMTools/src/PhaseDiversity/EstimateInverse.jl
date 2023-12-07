module EstimateInverse
using ..LatticeTools: Lattice, dualShiftLattice
using FFTW: fft, ifft, fftshift, ifftshift

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



end # module EstimateInverse