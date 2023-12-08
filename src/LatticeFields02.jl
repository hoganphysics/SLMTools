# LatticeFields.jl




#-------------------- Dual lattices ------------------------------

function dualLattice(L::Lattice{N},flambda::Number=1) where N
    # The step size of the dual lattice is flambda/(length of lattice) = flambda/(N * dx).
    # The periodicity (~largest element) of the dual lattice is twice the Nyquist frequency, ~flambda/dx
    return ( ( (0:length(l)-1) .* flambda ./ (length(l)*step(l)) for l in L )..., )
end

function dualShiftLattice(L::Lattice{N},flambda::Number=1) where N
    # These frequencies are equivalent to those of fftshift(dualLattice), but the initial entries are negative
    return ( ( (-floor(Int,length(l)/2):floor(Int,(length(l)-1)/2)) .* flambda ./ (length(l)*step(l)) for l in L )..., )
end





#---------------------- Fourier helpers --------------------------------------

# Lattice duality query
ldq(L1::Lattice,L2::Lattice,flambda=1) = (all(isapprox.(dualShiftLattice(L1,flambda),L2)) || throw(DomainError((x.L,y.L), "Unequal lattices.")); return nothing;)
ldq(f1::LF,f2::LF) = (all(isapprox.(dualShiftLattice(f1.L,f1.flambda),f2.L)) && f1.flambda == f2.flambda || throw(DomainError((f1.L,f2.L, f1.flambda, f2.flambda), "Non-dual lattices.")); return nothing;)

function latticeDisplacement(L::Lattice)
    return [l[1] + floor(length(l)/2)*step(l) for l in L]
end

function toDim(v,d::Int,n::Int)
    reshape(v,((i==d ? length(v) : 1) for i=1:n)...)
end

function dualPhase(L::Lattice{N}; dL = nothing) where N
    if isnothing(dL)
        dL = dualShiftLattice(L)
    end
    b = latticeDisplacement(L)
    return LF{ComplexPhase}(exp.(-2*pi*im .* .+((b[i] * toDim(dL[i],i,N) for i=1:N)...)) , dL)
end

