module DualLattices
using ..LatticeCore
using Interpolations: LinearInterpolation  
linear_interpolation = LinearInterpolation
using FFTW: fft, ifft, fftshift, ifftshift

export dualLattice, dualShiftLattice, dualPhase, dualate, ldq

export Lattice, elq, RealPhase, Generic, Phase, FieldVal, ComplexPhase, UPhase, UnwrappedPhase, S1Phase
export Intensity, Amplitude, Modulus, RealAmplitude, RealAmp, ComplexAmplitude, ComplexAmp, LatticeField, LF
export subfield, wrap, square, sublattice

# export natlat, padout, sft, isft, latticeDisplacement, toDim


"""
    dualLattice(L::Lattice{N}, flambda::Number=1) where N

    Generally this is not the function you want--you probably want "dualShiftLattice".  This function creates the dual lattice of a given lattice `L` with all positive frequencies.  "dualShiftLattice" creates the corresponding dual lattice with the latter half of the frequencies aliased to negative values, which is the form in which the FFT approximates the continuous Fourier transform. The parameter `flambda` represents the product of the wavelength and the lens focal length used in the Fourier transform. The step size of the dual lattice is `flambda / (N * dx)`, with `N` being the size of the original lattice and `dx` its step size. The periodicity of the dual lattice is approximately twice the Nyquist frequency.

    Arguments:
    - L::Lattice{N}: The original lattice.
    - flambda::Number: Product of the wavelength and lens focal length (default is 1).

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

    Create a centered dual lattice of a given lattice `L`. The `flambda` parameter is the product of the wavelength and the lens focal length used in the Fourier transform. This function generates a lattice akin to the `fftshift` of dualLattice, with initial entries being negative.

    Arguments:
    - L::Lattice{N}: The original lattice.
    - flambda::Number: A product of the wavelength and lens factor (default is 1).

    Returns:
    - Tuple: A tuple representing the centered dual lattice.
    """
function dualShiftLattice(L::Lattice{N}, flambda::Number=1) where {N}
    # These frequencies are equivalent to those of fftshift(dualLattice), but the initial entries are negative
    return (((-floor(Int, length(l) / 2):floor(Int, (length(l) - 1) / 2)) .* flambda ./ (length(l) * step(l)) for l in L)...,)
end
LatticeCore.sft(v::LF{ComplexAmplitude}) = LF{ComplexAmplitude}(fftshift(fft(ifftshift(v.data))), dualShiftLattice(v.L, v.flambda), v.flambda)
LatticeCore.isft(v::LF{ComplexAmplitude}) = LF{ComplexAmplitude}(fftshift(ifft(ifftshift(v.data))), dualShiftLattice(v.L, v.flambda), v.flambda)

"""
    ldq(L1::Lattice, L2::Lattice, flambda=1)

    Check if two lattices, `L1` and `L2`, are duals of each other, considering a wavelength scaling factor `flambda`. The function computes the dual of `L1` using dualShiftLattice and compares it with `L2`.

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

    Check if the lattices of two `LatticeField` instances, `f1` and `f2`, are dual to each other. Also checks that they have the same `flambda` value.

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
    dualPhase(L::Lattice{N}, flambda::Real=1; dL = nothing) where N

    Calculate the dual phase of a given lattice `L` using a specified wavelength scaling factor `flambda`.  Uses supplied dual lattice dL if provided, otherwise computes it from dualShiftLattice. Returns a LatticeField{RealPhase}.

    Arguments:
    - `L::Lattice{N}`: The original lattice.
    - `flambda::Real`: Scale factor.

    Returns:
    - LatticeField{RealPhase,FloatF64, N}: LatticeField{RealPhase} representing the dual phase of the lattice `L`.

    The dual phase is computed using the formula `.+((b[i] * toDim(dL[i],i,N) for i=1:N)...) ./ flambda`, where `dL` is the shifted dual lattice and `b` is the lattice displacement.
    """
function dualPhase(L::Lattice{N}, flambda::Real=1.0; dL = nothing) where N
    if isnothing(dL)
        dL = dualShiftLattice(L)
    end
    b = latticeDisplacement(L)
    return LF{RealPhase}(.+((b[i] * toDim(dL[i],i,N) for i=1:N)...) ./ flambda , dL, flambda)
end


end # module DualLattices