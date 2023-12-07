



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
function getCostMatrix(L::Lattice{N};normalization=maximum) where N
    M = [sum( (L[k][I[k]] - L[k][J[k]])^2 for k=1:N ) for I in CartesianIndices(length.(L))[:], J in CartesianIndices(length.(L))[:]]
    return M ./ normalization(M)
end

"""
    getCostMatrix(Lμ::Lattice{N}, Lν::Lattice{N}; normalization=maximum) where N

    Calculate a cost matrix between two lattices `Lμ` and `Lν`.

    # Arguments
    - `Lμ::Lattice{N}`, `Lν::Lattice{N}`: Two lattices represented as multi-dimensional arrays.
    - `normalization`: A function used for normalizing the cost matrix.

    # Returns
    - A normalized cost matrix representing the cost between points in `Lμ` and `Lν`.

    This function computes the cost matrix based on the Euclidean distance squared between corresponding points in the two lattices.
    """
function getCostMatrix(Lμ::Lattice{N},Lν::Lattice{N}; normalization=maximum) where N
    M = [sum( (Lμ[k][I[k]] - Lν[k][J[k]])^2 for k=1:N ) for I in CartesianIndices(length.(Lμ))[:], J in CartesianIndices(length.(Lν))[:]]
    return M ./ normalization(M)
end

"""
    pdCostMatrix(L::Lattice{N}, αRoot::Real, αTarget::Real; normalization=maximum) where {N}

    Calculate the cost matrix for optimal transport problem.

    # Arguments
    - `L::Lattice{N}`: The lattice on which the cost is calculated.
    - `αRoot::Real`: The root alpha value.
    - `αTarget::Real`: The target alpha value.
    - `normalization`: The normalization function to apply to the cost matrix. Default is `maximum`.

    # Returns
    - `M`: The normalized cost matrix.
    """
function pdCostMatrix(L::Lattice{N}, αRoot::Real, αTarget::Real; normalization=maximum) where {N}
    Δ = αRoot - αTarget
    M = [sum((L[k][I[k]] - L[k][J[k]] / Δ)^2 for k = 1:N) for I in CartesianIndices(length.(L))[:], J in CartesianIndices(length.(L))[:]]
    return M ./ normalization(M)
end

"""
    pdotBeamEstimate(G2Root::AbstractArray{<:Number,N}, G2Target::AbstractArray{<:Number,N}, αRoot::Real, αTarget::Real, L::Lattice{N}, ε::Real; options...) where {N}

    Estimate the beam using optimal transport.

    # Arguments
    - `G2Root::AbstractArray{<:Number,N}`: The root distribution.
    - `G2Target::AbstractArray{<:Number,N}`: The target distribution.
    - `αRoot::Real`: The root alpha value.
    - `αTarget::Real`: The target alpha value.
    - `L::Lattice{N}`: The lattice on which the cost is calculated.
    - `ε::Real`: The regularization parameter for the Sinkhorn algorithm.
    - `options...`: Additional options for the Sinkhorn algorithm.

    # Returns
    - The estimated beam.

    # Errors
    - If the distributions do not have the same shape, an error is thrown.
    - If the Sinkhorn algorithm returns NaN, an error is thrown.
    """
function pdotBeamEstimate(G2Root::AbstractArray{<:Number,N}, G2Target::AbstractArray{<:Number,N}, αRoot::Real, αTarget::Real, L::Lattice{N}, ε::Real; options...) where {N}
    size(G2Root) == size(G2Target) == length.(L) || error("Distributions must have same shape.")
    C = pdCostMatrix(L, αRoot, αTarget)
    γ = sinkhorn(G2Root[:], G2Target[:], C, ε; options...)
    any(isnan, γ) && error("sinkhorn returned nan.  Try reducing epsilon.")
    Γ = mapify(γ, (length.(L)..., length.(L)...), L) ./ (αRoot - αTarget)
    Φ = scalarPotentialN(Γ, L)
    CI = CartesianIndices(length.(L))
    Φ .-= [sum(L[i][I[i]]^2 for i = 1:N) / (2 * (αRoot - αTarget)) for I in CI]
    dL = dualShiftLattice(L, 1)
    return isft(sqrt.(G2Root) .* exp.(2π * im * Φ)) .* [exp(-2π * im * αRoot / 2 * sum(L[i][I[i]]^2 for i = 1:N)) for I in CI]
end


"""
    safeInverse(x::Number)

    Calculate the multiplicative inverse of `x` safely.

    # Arguments
    - `x::Number`: The number to invert.

    # Returns
    - The multiplicative inverse of `x` if `x` is non-zero, otherwise returns zero of the same type as `x`.

    This function is a safeguard against division by zero errors by checking if `x` is zero before attempting to calculate its inverse.
    """
safeInverse(x::Number) = iszero(x) ? zero(typeof(x)) : 1/x

"""
    mapify(plan::AbstractArray{<:Number,2}, Lμ::Lattice{N}, Lν::Lattice{N}) where N

    Transform a transportation plan into a vector field.

    # Arguments
    - `plan::AbstractArray{<:Number,2}`: A matrix representing the transportation plan between two distributions.
    - `Lμ::Lattice{N}`: The source lattice.
    - `Lν::Lattice{N}`: The target lattice.

    # Returns
    - A vector field representing the transportation from `Lμ` to `Lν` as prescribed by `plan`.

    This function reshapes the transportation plan into a higher-dimensional array to align with the lattices. It then computes the vector field by accumulating the plan over the target lattice and normalizing it by the sum over the corresponding dimensions.
    """
function mapify(plan::AbstractArray{<:Number,2},Lμ::Lattice{N},Lν::Lattice{N}) where N
    sμ, sν = length.(Lμ), length.(Lν)
    lμ, lν = size(plan)
    p = reshape(plan,(lμ,sν...))
    vf = zeros(lμ,N)    # Vector field output.  Last coordinate indexes the vector components. 
    for i=1:N
        x = reshape(sum(p,dims=[(2:i)...,(i+2:N+1)...]),(lμ,sν[i]))
        vf[:,i] = x * Lν[i] .* safeInverse.(sum(x,dims=2))
    end
    return reshape(vf,(sμ...,N))
end

"""
    scalarPotential2(A::AbstractArray{T,3}, L::Lattice{2}; idx=nothing) where {T<:Number}

    Compute the scalar potential for a 2D vector field.

    # Arguments
    - `A::AbstractArray{T,3}`: A 3D array representing a 2D vector field, with the last dimension indexing the vector components.
    - `L::Lattice{2}`: A 2D lattice representing the spatial domain of the vector field.
    - `idx`: An optional Cartesian index specifying the reference point for the potential calculation.

    # Returns
    - A 2D array representing the scalar potential of the vector field `A` over the lattice `L`.

    This function computes a discrete scalar potential of a 2D vector field by integrating the field along each axis. The potential is set to zero at the `idx` point, or at the center of the lattice if `idx` is not provided. The function currently only supports 2D fields.
    """
function scalarPotential2(A::AbstractArray{T,3},L::Lattice{2};idx=nothing) where {T<:Number}
    # This only works in 2D currently. 
    #    M == N+1 || error("Dimension mismatch.")
    CI = CartesianIndices(length.(L))
    isnothing(idx) && (idx = CartesianIndex((length.(L).÷2)...))
    y = cumsum(A[idx[1]:idx[1],:,2],dims=2) * step(L[2])
    y .-= y[idx[2]]
    x = cumsum(A[:,:,1],dims=1) * step(L[1])
    x .+= y - x[idx[1]:idx[1],:]
    return x
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
function scalarPotentialN(A::AbstractArray{T,M},L::Lattice{N}; idx = nothing, dimOrder=1:N) where {M,T<:Number,N}
    M == N+1 && N == size(A)[end] || error("A must have dimension one greater than L, and size(A)[end] must equal the dimension of L.")
    isnothing(idx) && (idx = CartesianIndex((length.(L).÷2)...))
    CI = CartesianIndices(size(A)[1:end-1])
    return .+( (hyperSum(A[CI,dimOrder[i]],idx,dimOrder[i],dimOrder[1:(i-1)])*step(L[dimOrder[i]]) for i=1:N)... )
end

"""
    hyperSum(A::AbstractArray{T,N}, originIdx::CartesianIndex, sumDim::Int, fixDims) where {T<:Number, N}

    Perform a cumulative sum of array `A` in a hyperplane fixed by `fixDims` at `originIdx`.

    # Arguments
    - `A::AbstractArray{T,N}`: An N-dimensional array to sum over.
    - `originIdx::CartesianIndex`: The index of the origin in the cumulative sum.
    - `sumDim::Int`: The dimension along which the cumulative sum is performed.
    - `fixDims`: The dimensions to be fixed (i.e., not summed over).

    # Returns
    - An array with the same dimensions as `A`, containing the cumulative sum over the specified dimension, adjusted by the value at `originIdx`.

    This function performs a cumulative sum over the dimension `sumDim` while fixing the other specified dimensions. It adjusts the sum relative to the value at `originIdx`, effectively setting the potential to zero at that point.
    """
function hyperSum(A::AbstractArray{T,N},originIdx::CartesianIndex,sumDim::Int,fixDims) where {T<:Number,N}
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

"""
    OTphase(μ::AbstractArray{<:Number,N}, ν::AbstractArray{<:Number,N}, ε::Real; options...) where N

    Calculate the optimal transport phase between two distributions.

    # Arguments
    - `μ::AbstractArray{<:Number,N}`: A distribution array.
    - `ν::AbstractArray{<:Number,N}`: Another distribution array of the same shape as `μ`.
    - `ε::Real`: A regularization parameter for the Sinkhorn algorithm.

    # Keyword Arguments
    - `options`: Additional options for the Sinkhorn algorithm.

    # Returns
    - The scalar potential array representing the optimal transport phase between `μ` and `ν`.

    This function first ensures that the input distributions have the same shape, then constructs a cost matrix based on a natural lattice. It uses the Sinkhorn algorithm to find the optimal transport plan, which is then used to compute the scalar potential.
    """
function OTphase(μ::AbstractArray{<:Number,N},ν::AbstractArray{<:Number,N},ε::Real;options...) where N
    
    size(μ) == size(ν) || error("Distributions must have same shape.")
    L = ( (natlat(n) for n in size(μ))..., )
    C = getCostMatrix(L)
    γ = sinkhorn(μ[:], ν[:], C, ε; options...)
    any(isnan,γ) && error("sinkhorn returned nan.  Try changing epsilon.")
    Γ = mapify(γ,L,L)
    Φ = scalarPotentialN(Γ,L)
    return Φ
end

"""
    OTphase(μ::AbstractArray{<:Number,N}, Lμ::Lattice{N}, ν::AbstractArray{<:Number,N}, Lν::Lattice{N}, ε::Real; options...) where N

    Calculate the optimal transport phase between two distributions on specified lattices.

    # Arguments
    - `μ::AbstractArray{<:Number,N}`: A distribution array corresponding to the `Lμ` lattice.
    - `Lμ::Lattice{N}`: The lattice corresponding to distribution `μ`.
    - `ν::AbstractArray{<:Number,N}`: Another distribution array corresponding to the `Lν` lattice.
    - `Lν::Lattice{N}`: The lattice corresponding to distribution `ν`.
    - `ε::Real`: A regularization parameter for the Sinkhorn algorithm.

    # Keyword Arguments
    - `options`: Additional options for the Sinkhorn algorithm.

    # Returns
    - The scalar potential array representing the optmial transport phase between the distributions `μ` and `ν`.

    The function checks for size consistency between the distributions and their corresponding lattices. It then calculates the cost matrix for the given lattices, applies the Sinkhorn algorithm to find the optimal transport plan, and computes the scalar potential.
    """
function OTphase(μ::AbstractArray{<:Number,N},Lμ::Lattice{N},ν::AbstractArray{<:Number,N},Lν::Lattice{N},ε::Real;options...) where N
    (size(μ) ==  length.(Lμ)) && (size(ν) == length.(Lν)) || error("Mismatched disribution and lattice sizes.")
    C = getCostMatrix(Lμ,Lν)
    γ = sinkhorn(μ[:],ν[:],C,ε; options...)
    any(isnan,γ) && error("sinkhorn returned nan.  Try changing epsilon.")
    Γ = mapify(γ,Lμ,Lν)
    Φ = scalarPotentialN(Γ, Lμ; dimOrder=N:-1:1)    # Formerly scalarPotential1.  Now will break in dimensions other than 2. 
    return Φ
end

"""
    OTreducedPhase(μ::AbstractArray{<:Number,N}, ν::AbstractArray{<:Number,N}, ε::Real, r::Int; options...) where N

    Compute the optical transport phase for downsampled distributions.

    # Arguments
    - `μ::AbstractArray{<:Number,N}`: A distribution array.
    - `ν::AbstractArray{<:Number,N}`: Another distribution array of the same shape as `μ`.
    - `ε::Real`: A regularization parameter for the Sinkhorn algorithm.
    - `r::Int`: Downsampling factor.

    # Keyword Arguments
    - `options`: Additional options for the Sinkhorn algorithm.

    # Returns
    - An array representing the interpolated optical transport phase for the original distributions.

    This function downsamples the input distributions by a factor of `r`, calculates the optical transport phase for the downsampled distributions, and then interpolates this phase back to the original resolution using cubic spline interpolation.
    """
function OTreducedPhase(μ::AbstractArray{<:Number,N},ν::AbstractArray{<:Number,N},ε::Real,r::Int;options...) where N
    size(μ) == size(ν) || error("Distributions must have same shape.")
    size(μ) .% r == 0 .* size(μ) || error("Downsample factor must divide distribution size.")
    U = downsample(μ,r)
    V = downsample(ν,r)
    ΦD = OTphase(U,V,ε;options...)
    L = ( (natlat(n) for n in size(μ))..., )
    LD = ( (l[1:r:end] .+ (l[r]-l[1])/2 for l in L)..., )
    Φ = cubic_spline_interpolation(LD,ΦD,extrapolation_bc=Flat())[L...] .* r
    return Φ
end

