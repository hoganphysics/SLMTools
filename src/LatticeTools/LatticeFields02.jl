# LatticeFields.jl




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

