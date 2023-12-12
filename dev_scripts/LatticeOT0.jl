# LatticeOT.jl
# Optimal transport methods for LatticeFields
# Requires: OptimalTransport.jl, FFTW
# There's just one place where FFTW is required, and this could reasonably be eliminated by moving that function elsewhere. 

#-------------- Cost matrices -------------------------------------------

function getCostMatrix(L::Lattice{N};normalization=maximum) where N
    M = [sum( (L[k][I[k]] - L[k][J[k]])^2 for k=1:N ) for I in CartesianIndices(length.(L))[:], J in CartesianIndices(length.(L))[:]]
    return M ./ normalization(M)
end

function getCostMatrix(Lμ::Lattice{N},Lv::Lattice{N}; normalization=maximum) where N
    M = [sum( (Lμ[k][I[k]] - Lv[k][J[k]])^2 for k=1:N ) for I in CartesianIndices(length.(Lμ))[:], J in CartesianIndices(length.(Lv))[:]]
    return M ./ normalization(M)
end

function pdCostMatrix(LRoot::Lattice{N},LTarget::Lattice{N},αRoot::Real,αTarget::Real; normalization=maximum) where N
    Δ = αRoot - αTarget
    M = [sum( (LRoot[k][I[k]] - LTarget[k][J[k]] / Δ)^2 for k=1:N ) for I in CartesianIndices(length.(LRoot))[:], J in CartesianIndices(length.(LTarget))[:]]
    return M ./ normalization(M)
end

#--------------------- Misc ----------------------------------------------

function getMoment(A::AbstractArray{<:Number,N},m::NTuple{M,Integer}) where {N,M}
    # Returns the moment <x_{i_1} x_{i_2} ... x_{i_M}> of A, where m = (i_1,i_2,...,i_M)
    pows = zeros(Int,N)
    for j in m
        pows[j]+=1
    end
    nondims = filter(x->pows[x]==0,1:N)
    indims = filter(x->pows[x]!=0,1:N)
    B = sum(A,dims=nondims)
    return sum( B .* [prod([I[j]^pows[j] for j in indims]) for I in CartesianIndices(size(B))] )
end;

#---------------------- Plan to map conversion -----------------------------------------

safeInverse(x::Number) = iszero(x) ? zero(typeof(x)) : 1/x

function mapify(plan::AbstractArray{<:Number,2},Lμ::Lattice{N},Lv::Lattice{N}) where N
    sμ, sv = length.(Lμ), length.(Lv)
    lμ, Lv = size(plan)
    p = reshape(plan,(lμ,sv...))
    vf = zeros(lμ,N)    # Vector field output.  Last coordinate indexes the vector components. 
    for i=1:N
        x = reshape(sum(p,dims=[(2:i)...,(i+2:N+1)...]),(lμ,sv[i]))
        vf[:,i] = x * Lv[i] .* safeInverse.(sum(x,dims=2))
    end
    return reshape(vf,(sμ...,N))
end

#------------------------- Gradient integration ----------------------------------------
# At some point, this should be replaced with a Helmholtz decomposition.  This is a known source of error, 
# especially for complex distributions, such as those typically output by phase diversity on real data. 

function hyperSum(A::AbstractArray{T,N},originIdx::CartesianIndex,sumDim::Int,fixDims) where {T<:Number,N}
    # Slices array A along dimensions fixDims and cumsums along dimension sumDim.  Does this for each value of 
    # the remaining dimensions.  originIdx indicates which CartesianIndex should be zero in the cumsum.
    sumDim in fixDims && error("sumDim can't be in fixDims")
    sumDim in 1:N && all(d in 1:N for d in fixDims) || error("More dimensions than array axes.")

    dims = (( (i in fixDims || i == sumDim ? size(A)[i] : (originIdx[i]:originIdx[i]))  for i=1:N)...,)
    CI1 = CartesianIndices(dims)
    sumBox = cumsum(A[CI1],dims=sumDim)
    
    dims = (( (i in fixDims ? size(A)[i] : (i==sumDim ? (originIdx[i]:originIdx[i]) : 1))  for i=1:N)...,)
    CI2 = CartesianIndices(dims)
    originSlice = sumBox[CI2]
    return sumBox .- originSlice
end

function scalarPotentialN(A::AbstractArray{T,M},L::Lattice{N}; idx = nothing, dimOrder=1:N) where {M,T<:Number,N}
    # Computes the scalar potential of vector field A, where the last dimension of A stores vector components. 
    # idx is the origin of the summation, i.e. where the potential is zero.  dimOrder is the order in which dimensions are 
    # summed, and unfortunately matters.  
    M == N+1 && N == size(A)[end] || error("A must have dimension one greater than L, and size(A)[end] must equal the dimension of L.")
    isnothing(idx) && (idx = CartesianIndex((length.(L).÷2)...))
    CI = CartesianIndices(size(A)[1:end-1])
    return .+( (hyperSum(A[CI,dimOrder[i]],idx,dimOrder[i],dimOrder[1:(i-1)])*step(L[dimOrder[i]]) for i=1:N)... )
end


#----------------------- OT solver -------------------------------------------------

function normalizeDistribution(U::AbstractArray{T,N}) where {T<:Number,N}
    # Makes an arbitrary array into a positive real array with sum 1, i.e. a probability distribution. 
    U = abs.(U)
    return U ./ sum(U)
end

function otPhase(U::LatticeField{Intensity,<:Real,N},V::LatticeField{Intensity,<:Number,N},ε::Real;options...) where N
    # Solves the optimal transport problem between distributions U and V, returning a LatticeField{RealPhase} representing
    # the scalar potential of the transport map. 
    u,v = normalizeDistribution(U.data), normalizeDistribution(V.data)
    C = getCostMatrix(U.L,V.L)
    γ = sinkhorn(u[:],v[:],C,ε; options...)
    any(isnan,γ) && error("sinkhorn returned nan.  Try changing epsilon.")
    Γ = mapify(γ,U.L,V.L)
    Φ = scalarPotentialN(Γ,U.L;dimOrder=N:-1:1) 
    return LF{RealPhase}(Φ,U.L)
end

#------------------------ Phase diversity OT solver --------------------------------

function pdotPhase(G2Root::LatticeField{Intensity,<:Real,N},G2Target::LatticeField{Intensity,<:Real,N},αRoot::Real,αTarget::Real, βRoot::Vector, βTarget::Vector, ε::Real;options...) where N
    u,v = normalizeDistribution(G2Root.data), normalizeDistribution(G2Target.data)
    C = pdCostMatrix(G2Root.L,G2Target.L,αRoot,αTarget)
    γ = sinkhorn(u[:], v[:], C, ε; options...)
    any(isnan,γ) && error("sinkhorn returned nan.  Try changing epsilon.")
    Γ = mapify(γ,G2Root.L,G2Target.L) ./ (αRoot-αTarget)
    Φ = scalarPotentialN(Γ,G2Root.L)
    CI = CartesianIndices(length.(G2Root.L))
    dβ = βRoot .- βTarget
    Φ .-= [ sum( (G2Root.L[i][I[i]] - dβ[i])^2 for i=1:N)/(2*(αRoot-αTarget)) for I in CI ]
    return LF{RealPhase}(Φ,G2Root.L)
end

function pdotBeamEstimate(G2Root::LatticeField{Intensity,<:Real,N},G2Target::LatticeField{Intensity,<:Real,N}, αRoot::Real,αTarget::Real, βRoot::Vector, βTarget::Vector, ε::Real;options...) where N
    Φ = pdotPhase(G2Root,G2Target,αRoot,αTarget,βRoot,βTarget,ε;options...)
    dL = dualShiftLattice(G2Root.L)
    beam = isft(sqrt.(G2Root.data) .* exp.(2π*im*Φ)) .* [exp(-2π*im*αRoot/2 * sum(dL[i][I[i]]^2 for i=1:N)) for I in CI]
    return LF{ComplexAmplitude}(beam,dL)
end
