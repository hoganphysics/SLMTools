module DualLattices
using ..LatticeCore: Lattice, sft, isft
using Interpolations: linear_interpolation 

export dualLattice, dualShiftLattice, dualate

"""
    dualLattice(L::Lattice{N}, flambda::Number=1) where N

    Create the dual lattice of a given lattice `L`. The parameter `flambda` represents the product of the wavelength and the lens used in the Fourier transform. The step size of the dual lattice is `flambda / (N * dx)`, with `N` being the size of the original lattice and `dx` its step size. The periodicity of the dual lattice is approximately twice the Nyquist frequency.

    Arguments:
    - L::Lattice{N}: The original lattice.
    - flambda::Number: A product of the wavelength and lens factor (default is 1).

    Returns:
    - Tuple: A tuple representing the dual lattice.
    """
function dualLattice(L::Lattice{N}, flambda::Number=1) where {N}
    # The step size of the dual lattice is flambda/(length of lattice) = flambda/(N * dx).
    # The periodicity (~largest element) of the dual lattice is twice the Nyquist frequency, ~flambda/dx
    return (((0:length(l)-1) .* flambda ./ (length(l) * step(l)) for l in L)...,)
end

"""
    dualShiftLattice(L::Lattice{N}, flambda::Number=1) where N

    Create a frequency-shifted dual lattice of a given lattice `L`. The `flambda` parameter is the product of the wavelength and the lens used in the Fourier transform. This function generates a lattice akin to the `fftshift` of the regular dual lattice, with initial entries being negative.

    Arguments:
    - L::Lattice{N}: The original lattice.
    - flambda::Number: A product of the wavelength and lens factor (default is 1).

    Returns:
    - Tuple: A tuple representing the frequency-shifted dual lattice.
    """
function dualShiftLattice(L::Lattice{N}, flambda::Number=1) where {N}
    # These frequencies are equivalent to those of fftshift(dualLattice), but the initial entries are negative
    return (((-floor(Int, length(l) / 2):floor(Int, (length(l) - 1) / 2)) .* flambda ./ (length(l) * step(l)) for l in L)...,)
end




""" dualate(imgs, slmGrid, camGrid, θ, flambda; roi=nothing)

    This function performs a dualate operation on the given images.

    # Arguments
    - `imgs`: An array of images to be processed.
    - `slmGrid`: The SLM grid for the images.
    - `camGrid`: The camera grid for the images. This should have physical units, not pixel units.
    - `θ`: The angle output of sweepCoords.
    - `flambda`: The scaling factor for the dual lattice.
    - `roi`: (optional) The region of interest in the images. If not provided, the entire image is used.

    # Returns
    - An array of images that have been processed with the dualate operation.

    # Errors
    - If the camera grid and image size are unequal, an error is thrown."""
function dualate(imgs, slmGrid, camGrid, θ, flambda; roi=nothing)
    if !isnothing(roi)
        if size(imgs[1]) != length.(roi)
            imgs = [i[roi...] for i in imgs]
        end
        if length.(camGrid) != length.(roi)
            camGrid = ((camGrid[i][roi[i]] for i = 1:length(camGrid))...,)
        end
    end
    size(imgs[1]) == length.(camGrid) || error("Camera grid and image size unequal.")
    interps = [linear_interpolation(camGrid, imgs[j], extrapolation_bc=zero(typeof(imgs[1][1]))) for j in 1:length(imgs)]
    dL = dualShiftLattice(slmGrid, flambda)
    Lx0, Lxs, Ly0, Lys = dL[1][1], step(dL[1]), dL[2][1], step(dL[2])
    R = [cos(θ) -sin(θ); sin(θ) cos(θ)]
    r0, dx, dy = (R * [Lx0, Ly0]), Lxs * R[:, 1], Lys * R[:, 2]

    return [[interp((r0 + dx * I[1] + dy * I[2])...) for I in CartesianIndices(length.(dL))] for interp in interps]
end






















end # module DualLattices