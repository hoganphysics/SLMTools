module DualLattices
using ..LatticeCore
using Interpolations: linear_interpolation 

export dualLattice, dualShiftLattice, dualPhase, dualate, ldq

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


"""
    ldq(L1::Lattice, L2::Lattice, flambda=1)

    Check if two lattices, `L1` and `L2`, are duals of each other, considering a wavelength scaling factor `flambda`. The function computes the dual of `L1` using a shifted lattice approach and compares it with `L2`.

    Arguments:
    - `L1::Lattice`: The first lattice.
    - `L2::Lattice`: The second lattice, to be compared with the dual of `L1`.
    - `flambda`: A scaling factor for the wavelength (default is 1).

    Returns:
    - Nothing, but throws `DomainError` if lattices are not duals.

    Usage:
    ldq(lattice1, lattice2, flambda)
    """
ldq(L1::Lattice, L2::Lattice, flambda=1) = (all(isapprox.(dualShiftLattice(L1, flambda), L2)) || throw(DomainError((L1, L2), "Unequal lattices.")); return nothing)

"""
    ldq(f1::LF, f2::LF)

    Check if two `LatticeField` instances, `f1` and `f2`, are duals of each other. It compares the dual of the lattice of `f1` with the lattice of `f2` and checks if their `flambda` values are the same.

    Arguments:
    - `f1::LF`: The first LatticeField.
    - `f2::LF`: The second LatticeField.

    Returns:
    - Nothing, but throws `DomainError` if LatticeFields are not duals.

    Usage:
    ldq(latticeField1, latticeField2)
    """
ldq(f1::LF, f2::LF) = (all(isapprox.(dualShiftLattice(f1.L, f1.flambda), f2.L)) && f1.flambda == f2.flambda || throw(DomainError((f1.L, f2.L, f1.flambda, f2.flambda), "Non-dual lattices.")); return nothing)

"""
    dualPhase(L::Lattice{N}; dL = nothing) where N

    Create a `LatticeField` representing the phase factors for a dual lattice. This function is used in Fourier optics and other applications involving wave propagation and interference. The phase factors are calculated based on the dual lattice of `L` (or a provided `dL`) and the lattice displacement.

    Arguments:
    - `L::Lattice{N}`: The original lattice.
    - `dL`: An optional dual lattice. If not provided, the dual lattice is computed using `dualShiftLattice(L)`.

    Returns:
    - `LatticeField{ComplexPhase}`: A `LatticeField` with complex phase values representing the phase factors of the dual lattice.

    Usage:
    # Without specifying dL
    dualPhaseField = dualPhase(myLattice)

    # With specifying dL
    dualLattice = dualShiftLattice(myLattice)
    dualPhaseField = dualPhase(myLattice, dL=dualLattice)
    """
function dualPhase(L::Lattice{N}; dL=nothing) where {N}
    if isnothing(dL)
        dL = dualShiftLattice(L)
    end
    b = latticeDisplacement(L)
    return LF{ComplexPhase}(exp.(-2 * pi * im .* .+((b[i] * toDim(dL[i], i, N) for i = 1:N)...)), dL)
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