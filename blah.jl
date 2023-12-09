function dualPhase(L::Lattice{N}, dL::Lattice{N}=nothing) where {N}
    if isnothing(dL)
        dL = dualShiftLattice(L)
    end
    b = latticeDisplacement(L)
    return LF{ComplexPhase}(exp.(-2 * pi * im .* .+((b[i] * toDim(dL[i], i, N) for i = 1:N)...)), dL)
end


function dualPhase(L::Lattice{N}, flambda::Real=1) where {N}
    dL = dualShiftLattice(L, flambda)
    b = latticeDisplacement(L)
    return exp.(-2 * pi * im .* .+((b[i] * toDim(dL[i], i, N) for i = 1:N)...) ./ flambda)
end