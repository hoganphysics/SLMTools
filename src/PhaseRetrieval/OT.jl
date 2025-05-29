module OT

using ..LatticeTools
using ..Misc
using OptimalTransport: sinkhorn
using Interpolations: Linear
using Images: Gray

export getCostMatrix, pdCostMatrix, pdotPhase, pdotBeamEstimate, mapify, scalarPotentialN, otPhase, hyperSum2
function look(f::AbstractArray{T,N}) where {T<:Real,N}
    return Gray.(f ./ maximum(f))
end
"""
    getCostMatrix(L::Lattice{N}; normalization=maximum) where N

    Calculate a cost matrix for a given lattice `L`.

    # Arguments
    - `L::Lattice{N}`: A lattice represented as a multi-dimensional array.
    - `normalization`: A function used for normalizing the cost matrix.

    # Returns
    - A normalized cost matrix derived from the lattice `L`.

    The cost matrix is calculated based on the Euclidean distance squared between points in the lattice `L`.
    """
function getCostMatrix(L::Lattice{N}; normalization=maximum) where {N}
    M = [sum((L[k][I[k]] - L[k][J[k]])^2 for k = 1:N) for I in CartesianIndices(length.(L))[:], J in CartesianIndices(length.(L))[:]]
    return M ./ normalization(M)
end

"""
    getCostMatrix(Lμ::Lattice{N}, Lv::Lattice{N}; normalization=maximum) where N

    Calculate a cost matrix between two lattices `Lμ` and `Lv`.

    # Arguments
    - `Lμ::Lattice{N}`, `Lv::Lattice{N}`: Two lattices represented as multi-dimensional arrays.
    - `normalization`: A function used for normalizing the cost matrix.

    # Returns
    - A normalized cost matrix representing the cost between points in `Lμ` and `Lv`.

    This function computes the cost matrix based on the Euclidean distance squared between corresponding points in the two lattices.
    """
function getCostMatrix(Lμ::Lattice{N}, Lv::Lattice{N}; normalization=maximum) where {N}
    M = [sum((Lμ[k][I[k]] - Lv[k][J[k]])^2 for k = 1:N) for I in CartesianIndices(length.(Lμ))[:], J in CartesianIndices(length.(Lv))[:]]
    return M ./ normalization(M)
end

"""
    pdCostMatrix(L::Lattice{N}, αRoot::Real, αTarget::Real; normalization=maximum) where {N}

    DEPRECATED: This function is not necessary, as using getCostMatrix results in the same final answer
                and better performance of the sinkhorn method.  We once thought this function necessary,
                but later understood from scaling properties of OT that it is not.  It may still be 
                occasionally useful for visualization of the cost matrix, so we leave it here. 
    
    Calculate the cost matrix for optimal transport problem. 

    # Arguments
    - `L::Lattice{N}`: The lattice on which the cost is calculated.
    - `αRoot::Real`: The root alpha value.
    - `αTarget::Real`: The target alpha value.
    - `normalization`: The normalization function to apply to the cost matrix. Default is `maximum`.

    # Returns
    - `M`: The normalized cost matrix.
    """
function pdCostMatrix(LRoot::Lattice{N}, LTarget::Lattice{N}, αRoot::Real, αTarget::Real; normalization=maximum, flambda::Real=1.0) where {N}
    Δ = (αRoot - αTarget)*flambda
    M = [sum((LRoot[k][I[k]] - LTarget[k][J[k]] / Δ)^2 for k = 1:N) for I in CartesianIndices(length.(LRoot))[:], J in CartesianIndices(length.(LTarget))[:]]
    return M ./ normalization(M)
end


"""
    mapify(plan::AbstractArray{<:Number,2}, Lμ::Lattice{N}, Lv::Lattice{N}) where N

    Transform a transportation plan into a vector field.

    # Arguments
    - `plan::AbstractArray{<:Number,2}`: A matrix representing the transportation plan between two distributions.
    - `Lμ::Lattice{N}`: The source lattice.
    - `Lv::Lattice{N}`: The target lattice.

    # Returns
    - A vector field representing the transportation from `Lμ` to `Lv` as prescribed by `plan`.

    This function reshapes the transportation plan into a higher-dimensional array to align with the lattices. 
        It then computes the vector field by accumulating the plan over the target lattice and normalizing it by the sum over the corresponding dimensions.
    """
function mapify(plan::AbstractArray{<:Number,2},Lμ::Lattice{N},Lv::Lattice{N}) where N
    sμ, sv = length.(Lμ), length.(Lv)
    lμ, lv = size(plan)
    p = reshape(plan,(lμ,sv...))
    vf = zeros(lμ,N)    # Vector field output.  Last coordinate indexes the vector components. 
    for i=1:N
        x = reshape(sum(p,dims=[(2:i)...,(i+2:N+1)...]),(lμ,sv[i]))
        vf[:,i] = x * Lv[i] .* safeInverse.(sum(x,dims=2))
    end
    return reshape(vf,(sμ...,N))
end

"""
    hyperSum2(A::AbstractArray{T,N}, originIdx::CartesianIndex, sumDim::Int, fixDims) where {T<:Number, N}

    Perform a cumulative sum of array `A` in a hyperplane fixed by `fixDims` at `originIdx`.  
    Differs from hyperSum by adding a trapezoid rule correction term which makes this a higher order integration method. 

    # Arguments
    - `A::AbstractArray{T,N}`: An N-dimensional array to sum over.
    - `originIdx::CartesianIndex`: The index of the origin in the cumulative sum.
    - `sumDim::Int`: The dimension along which the cumulative sum is performed.
    - `fixDims`: The dimensions to be fixed (i.e., not summed over).

    # Returns
    - An array with the same dimensions as `A`, containing the cumulative sum over the specified dimension, adjusted by the value at `originIdx`.

    This function performs a cumulative sum over the dimension `sumDim` while fixing the other specified dimensions. It adjusts the sum relative to the value at `originIdx`, effectively setting the potential to zero at that point.
    """
function hyperSum2(A::AbstractArray{T,N}, originIdx::CartesianIndex, sumDim::Int, fixDims) where {T<:Number,N}
    # Slices array A along dimensions fixDims and cumsums along dimension sumDim.  originIdx indicates which CartesianIndex
    # should be zero in the cumsum. Collapses remaining dimensions to singletons at the respective originIdx coordinate.  
    # This function differs from hyperSum in that it adds a trapezoid rule correction term. This has a dramatic effect on the 
    # output of otPhase. 
    sumDim in fixDims && error("sumDim can't be in fixDims")
    sumDim in 1:N && all(d in 1:N for d in fixDims) || error("More dimensions than array axes.")

    dims = Tuple((i in fixDims || i == sumDim ? size(A)[i] : (originIdx[i]:originIdx[i])) for i = 1:N)
    CI1 = CartesianIndices(dims)
    sumBox = cumsum(A[CI1], dims=sumDim) .- A[CI1] ./ 2

    dims = Tuple((i in fixDims ? size(A)[i] : (i == sumDim ? (originIdx[i]:originIdx[i]) : 1)) for i = 1:N)
    CI2 = CartesianIndices(dims)
    originSlice = sumBox[CI2]
    return sumBox .- originSlice
end

"""
    scalarPotentialN(A::AbstractArray{T,M}, L::Lattice{N}; idx = nothing, dimOrder=1:N) where {M, T<:Number, N}

    Calculate the scalar potential for an N-dimensional vector field.

    # Arguments
    - `A::AbstractArray{T,M}`: A multi-dimensional array representing an N-dimensional vector field with `M = N+1`.
    - `L::Lattice{N}`: An N-dimensional lattice representing the spatial domain.
    - `idx`: An optional Cartesian index specifying the reference point for potential calculation.
    - `dimOrder`: An optional ordering of dimensions to be used in calculating the potential.

    # Returns
    - An N-dimensional array representing the scalar potential of the vector field `A`.

    This function verifies dimensionality consistency, then computes the scalar potential by summing over each dimension in the specified `dimOrder`, multiplying by the lattice spacing, and adjusting based on the reference point `idx`.
    """
function scalarPotentialN(A::AbstractArray{T,M}, L::Lattice{N}; idx=nothing, dimOrder=1:N) where {M,T<:Number,N}
    # Computes the scalar potential of vector field A, where the last dimension of A stores vector components. 
    # idx is the origin of the summation, i.e. where the potential is zero.  dimOrder is the order in which dimensions are 
    # summed, and unfortunately matters.  
    M == N + 1 && N == size(A)[end] || error("A must have dimension one greater than L, and size(A)[end] must equal the dimension of L.")
    isnothing(idx) && (idx = CartesianIndex((length.(L) .÷ 2)...))
    CI = CartesianIndices(size(A)[1:end-1])
    return .+((hyperSum2(A[CI, dimOrder[i]], idx, dimOrder[i], dimOrder[1:(i-1)]) * step(L[dimOrder[i]]) for i = 1:N)...)
end


function normalizeDistribution(U::AbstractArray{T,N}) where {T<:Number,N}
    # Makes an arbitrary array into a positive real array with sum 1, i.e. a probability distribution. 
    U = abs.(U)
    return U ./ sum(U)
end

"""
    otPhase(U::LatticeField{Intensity,<:Real,N}, V::LatticeField{Intensity,<:Number,N}, ε::Real; options...) where {N}

    Solve the optimal transport problem between two distributions represented by lattice fields.

    Given two distributions `U` and `V` represented as `LatticeField{Intensity}` objects, this function computes 
    the optimal transport map between them using the Sinkhorn algorithm. It returns a `LatticeField{RealPhase}` 
    representing the scalar potential of the transport map.

    # Arguments
    - `U::LatticeField{Intensity,<:Real,N}`: A lattice field representing the source distribution.
    - `V::LatticeField{Intensity,<:Number,N}`: A lattice field representing the target distribution.
    - `ε::Real`: The regularization parameter for the Sinkhorn algorithm.
    - `options...`: Additional options passed to the Sinkhorn algorithm.

    # Returns
    - `LatticeField{RealPhase}`: A lattice field representing the scalar potential of the optimal transport map.

    # Notes
    - The distributions `U` and `V` are first normalized to ensure they sum to 1.
    - The cost matrix for transport is computed based on the lattice parameters of `U` and `V`.
    - The Sinkhorn algorithm is used to compute the optimal transport plan, which is then converted into a scalar potential.

    # Errors
    - If `sinkhorn` returns `NaN`, an error is thrown suggesting to try changing the value of `ε`.

    # See Also
    - `normalizeDistribution`, `sinkhorn`, `scalarPotentialN`
    """
function otPhase(U::LatticeField{Intensity,<:Real,N},V::LatticeField{Intensity,<:Number,N},ε::Real;options...) where N
    # Solves the optimal transport problem between distributions U and V, returning a LatticeField{RealPhase} representing
    # the scalar potential of the transport map. 
	U.flambda == V.flambda || error("Unequal flambdas.")
    u,v = normalizeDistribution(U.data), normalizeDistribution(V.data)
    L1 = natlat(size(U))
    L2 = natlat(size(V))
    L2 = Tuple(L2[i] .* (maximum(L1[i])/maximum(L2[i])) for i=1:N)
    C = getCostMatrix(L1, L2)
    γ = sinkhorn(u[:],v[:],C,ε; options...)
    any(isnan,γ) && error("sinkhorn returned nan.  Try changing epsilon.")
    Γ = mapify(γ,U.L,V.L)
    Φ = scalarPotentialN(Γ,U.L;dimOrder=N:-1:1) / U.flambda
    return LF{RealPhase}(Φ,U.L,U.flambda)
end


"""
    pdotPhase(G2Root::LatticeField{Intensity,<:Real,N},
              G2Target::LatticeField{Intensity,<:Real,N},
              αRoot::Real,
              αTarget::Real,
              βRoot::Vector,
              βTarget::Vector,
              ε::Real; options...) where N
	
    # Returns
    - `LatticeField{RealPhase}`: A lattice field representing the phase associated with the LatticeField G2Root.

    # Errors
    - An error is thrown if `sinkhorn` returns `NaN`, suggesting a change in the value of `ε`.
    """

function pdotPhase(G2Root::LatticeField{Intensity,<:Real,N},G2Target::LatticeField{Intensity,<:Real,N},αRoot::Real,αTarget::Real, βRoot::Union{Vector,NTuple{N,Real}}, βTarget::Union{Vector,NTuple{N,Real}}, ε::Real;options...) where N
	G2Root.flambda == G2Target.flambda || error("Unequal flambdas.")
    u,v = normalizeDistribution(G2Root.data), normalizeDistribution(G2Target.data)
    L1 = natlat(size(G2Root))
    L2 = natlat(size(G2Target))
    L2 = Tuple(L2[i] .* (maximum(L1[i])/maximum(L2[i])) for i=1:N)
    C = getCostMatrix(L1, L2)
    γ = sinkhorn(u[:], v[:], C, ε; options...)
    any(isnan,γ) && error("sinkhorn returned nan.  Try changing epsilon.")
    Γ = mapify(γ,G2Root.L,G2Target.L) ./ (αRoot-αTarget)
    Φ = scalarPotentialN(Γ,G2Root.L)
    CI = CartesianIndices(length.(G2Root.L))
    dβ = βRoot .- βTarget
    Φ .-= [ sum( (G2Root.L[i][I[i]] - dβ[i])^2 for i=1:N)/(2*(αRoot-αTarget)) for I in CI ]
	Φ ./= G2Root.flambda^2
    return LF{RealPhase}(Φ,G2Root.L,G2Root.flambda)
end

"""
    pdotBeamEstimate(G2Root::LatticeField{Intensity,<:Real,N},G2Target::LatticeField{Intensity,<:Real,N}, αRoot::Real,αTarget::Real, βRoot::Vector, βTarget::Vector, ε::Real; LFine::Union{Nothing,Lattice{N}}=nothing, options...) where N
    Φ = pdotPhase(G2Root,G2Target,αRoot,αTarget,βRoot,βTarget,ε;options...) where N

    # Arguments
    - `G2Root::LatticeField{Intensity,<:Real,N}`: The larger of two diversity images
    - `G2Target::LatticeField{Intensity,<:Real,N}`: The smaller of two diversity images
    - `αRoot::Real`: The curvature (from file name ) of the larger diversity image
    - `αTarget::Real`: The curvature (from file name ) of the smaller diversity image
    - `βRoot::Vector`: The offset of the larger diversity image (linear phase)
    - `βTarget::Vector`: The offset of the smaller diversity image (linear phase)
    - `ε::Real`: The regularization parameter for the Sinkhorn algorithm.
    - `LFine::Union{Nothing,Lattice{N}}=nothing`: An optional lattice to upsample the result to (in the camera plane).

    # Returns
    - `LatticeField{ComplexAmplitude}`: A lattice field representing the inferred beam.
    """
function pdotBeamEstimate(G2Root::LatticeField{Intensity,<:Real,N},G2Target::LatticeField{Intensity,<:Real,N}, αRoot::Real,αTarget::Real, βRoot::Union{Vector,NTuple{N,Real}}, βTarget::Union{Vector,NTuple{N,Real}}, ε::Real; LFine::Union{Nothing,Lattice{N}}=nothing, options...) where N
    Φ = pdotPhase(G2Root,G2Target,αRoot,αTarget,βRoot,βTarget,ε;options...)
    GR = sqrt(G2Root)
    if !isnothing(LFine)
        GR = upsample(GR,LFine; bc=0) * upsample(Φ,LFine; bc=Linear())
    else
        GR = GR * Φ
    end
    dL = dualShiftLattice(GR.L,GR.flambda)
    CI = CartesianIndices(size(GR))
    divPhase = LF{RealPhase}([(αRoot/2 * sum(dL[i][I[i]]^2 for i=1:N) + sum(βRoot[i] * dL[i][I[i]] for i=1:N)) for I in CI],dL,GR.flambda)
    return isft(GR) * conj(divPhase)
end


# NEW CODE HERE:
function SinkhornConv2D(μ::Matrix{T}, ν::Matrix{T}, ϵ::Real, max_iter::Integer) where {T<:Real}
    
    # parameters
    N = size(μ)[1]
    u = ones(N, N)./N^2
    v = ones(N, N)./N^2 # this initialization doesn't matter

    # loss tracking
    loss = []
    prev = copy(u) 
    
    # create convolution matrix
    X, Y = natlat((N,N))
    A = exp.(-(X.^2 .+ Y'.^2) / (2 * N * ϵ))
    FA = sft(A)

    for i in range(1, max_iter)

        # row constraint
        row_sum = real.(isft(sft(u).*FA))
        v = ν ./ row_sum
        v = v ./ sum(v)
        
        # col constraint
        col_sum = real.(isft(sft(v).*FA))
        u = μ ./ col_sum
        u = u ./ sum(u)

        # log loss
        push!(loss, norm(u-prev))
        prev = copy(u)

    end
    
    # row constraint
    row_sum = real.(isft(sft(u).*FA))
    v = ν ./ row_sum
    
    # col constraint
    col_sum = real.(isft(sft(v).*FA))
    u = μ ./ col_sum
    
    return u, v, loss
end


function dualToGradients(u::Matrix{T}, v::Matrix{T}, μ::Matrix{S}, ϵ::Real) where {T<:Real, S<:Real}
    
    # create convolution matrix
    N = size(u)[1]
    X, Y = natlat(N,N)
    A = exp.(-(X.^2 .+ Y'.^2) / (2 * N * ϵ))
    FA = sft(A)
    
    # moment calculation
    XV = X .* v
    YV = Y' .* v

    # convolution step
    AXV = real.(isft(sft(XV).*FA))
    AYV = real.(isft(sft(YV).*FA))
    
    # compute gradients
    gradϕx = u .* AXV ./ μ
    gradϕy = u .* AYV ./ μ
    
    return gradϕx, gradϕy
end

function otQuickPhase(g2::LF{Intensity,T,N}, G2::LF{Intensity,T,N}, ϵ::Real, max_iter::Integer) where {T<:Real,N}
    u, v, loss = SinkhornConv2D(g2.data, G2.data, ϵ, max_iter) 
    gradϕx, gradϕy = dualToGradients(u, v, g2.data,  ϵ)
    Φ = scalarPotentialN(cat(gradϕx,gradϕy; dims=3), g2.L)
    return Φ
end 

end # module OTHelpers

