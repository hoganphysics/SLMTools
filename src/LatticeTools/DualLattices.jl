module DualLattices
using ..LatticeCore: Lattice, sft, isft

export dualLattice, dualShiftLattice, dualate


"""
    dualLattice(L::Lattice{N}, flambda::Number)
    Compute the dual lattice for a given lattice L in N dimensions and a scaling factor flambda.
    The dual lattice is used in frequency space and is related to the original lattice (L)
    by the Fourier transform. The step size of the dual lattice is determined by flambda
    and the spacing in the original lattice.
    Parameters:
    - L: The lattice for which to compute the dual.
    - flambda: The scaling factor for the dual lattice.
"""
function dualLattice(L::Lattice{N}, flambda::Number) where {N}
    return (((0:length(l)-1) .* flambda ./ (length(l) * step(l)) for l in L)...,)
end

"""
    dualShiftLattice(L::Lattice{N}, flambda::Number) 
    Compute a shifted version of the dual lattice.
    This function is similar to dualLattice but shifts the initial entries to be negative.
    This shift is often used in FFTs to center the zero-frequency component,
    making it analogous to the fftshift operation in frequency space."""
function dualShiftLattice(L::Lattice{N}, flambda::Number) where {N}
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