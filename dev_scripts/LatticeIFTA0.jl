# LatticeIFTA.jl
# Lattice versions of iterative Fourier transform algorithms for phase retrieval/generation and phase diversity/beam estimation.
# Requires: FFTW, LatticeFields.jl

#------------------ Gerchberg-Saxton ----------------------------------------------------------------
function gs(U::LF{Modulus,<:Real,N},V::LF{Modulus,<:Real,N},nit::Integer, Φ0::LF{ComplexPhase,<:Complex,N}) where N
    # Runs nit number of iterations of Gerchberg-Saxton on the distributions specified by U and V.  Φ0 is a starting phase guess.
    ldq(U,V), elq(U,Φ0)                              # Check for compatibility of the lattices
    ft = plan_fft(zeros(ComplexF64,size(U)))
    ift = plan_ifft(zeros(ComplexF64,size(V)))

    guess = ifftshift(U.data .* phasor.(Φ0.data))
    ushift = ifftshift(U.data)
    vshift = ifftshift(V.data)
    for i=1:nit
        guess = gsIter(guess,ushift,vshift,ft,ift)
    end
    return LF{ComplexPhase}(fftshift(phasor.(guess)),U.L)
end
phasor(z::ComplexF64) = iszero(z) ? one(ComplexF64) : z/abs(z)
function gsIter(guess::Array{Complex{Float64},N},u::Array{Float64, N},v::Array{Float64,N},ft::FFTW.cFFTWPlan,ift::AbstractFFTs.ScaledPlan) where N
    # Single iteration of Gerchberg-Saxton using provided Fourier operators
    return u .* phasor.(ift*(phasor.(ft*guess) .* v))
end

# Other methods.  These should just do the right thing on other input types. 
function gs(U::LF{Modulus,<:Real,N},V::LF{Modulus,<:Real,N},nit::Integer, Φ0::LF{RealPhase,<:Real,N}) where N
    return gs(U,V,nit,wrap(Φ0))
end
function gs(U::LF{Modulus,<:Real,N},V::LF{Modulus,<:Real,N},nit::Integer) where N
    return gs(U,V,nit,LF{RealPhase}(rand(size(U)...),U.L,U.flambda))
end
function gs(U::LF{Intensity,<:Real,N},V::LF{Intensity,<:Real,N},nit::Integer,Φ0::Union{Nothing,LF{ComplexPhase,<:Complex,N},LF{RealPhase,<:Real,N}}=nothing) where N
    return gs(sqrt(U),sqrt(V),nit,Φ0)
end


#--------------------- Phase diversity ----------------------------------------------

centroid(img::Array{T,N}) where {T,N} = [ sum( (1:size(img)[j]) .* sum(img,dims=delete!(Set(1:N),j))[:] ) for j=1:N ] ./ sum(img)
function window(img::Array{T,N}, w) where {T,N}
    c = CartesianIndex(round.(Int,centroid(img))...)
    return CartesianIndices(((1:w for i=1:N)...,)) .- CartesianIndex((floor(Int,w/2) for i=1:N)...) .+ c
end
function latticeDisplacement(L::Lattice)
    return [l[1] + floor(length(l)/2)*step(l) for l in L]
end
function toDim(v,d::Int,n::Int)
    reshape(v,((i==d ? length(v) : 1) for i=1:n)...)
end
function dualPhase(L::Lattice{N},flambda::Real=1) where N
    dL = dualShiftLattice(L,flambda)
    b = latticeDisplacement(L)
    return exp.(-2*pi*im .* .+((b[i] * toDim(dL[i],i,N) for i=1:N)...) ./ flambda)
end

function pdgs(imgs::NTuple{M,LF{Modulus,<:Number,N}},divPhases::NTuple{M,LF{<:Phase,<:Number,N}},nit::Integer,beamGuess::LF{ComplexAmp,<:Complex,N}) where {M,N}
    s = size(imgs[1])
    all(size(i)==s for i in imgs) && all(size(i)==size(beamGuess) for i in divPhases) || error("Input size mismatch.")
    [elq(d,beamGuess) for d in divPhases], [ldq(i,beamGuess) for i in imgs]
    ft = plan_fft(zeros(s))
    ift = plan_ifft(zeros(s))
    guess = ifftshift(beamGuess.data)
    phis = ((ifftshift(wrap(divPhases[i]).data .* dualPhase(imgs[i].L,imgs[i].flambda))*1.0 for i=1:M)...,)::NTuple{M,Array{ComplexF64,N}}
    mods = ((ifftshift(i.data)*1.0 for i in imgs)...,)::NTuple{M,Array{Float64,N}}
    for j=1:nit
        guess = pdgsIter(guess,phis,mods,ft,ift)
    end
    return LF{ComplexAmp}(fftshift(guess),beamGuess.L,beamGuess.flambda)
end
function pdgsIter(guess::Array{Complex{Float64},N},phis::NTuple{M,Array{ComplexF64,N}},mods::NTuple{M,Array{Float64,N}},
        ft::FFTW.cFFTWPlan,ift::AbstractFFTs.ScaledPlan) where {M,N}
    # Single iteration of Gerchberg-Saxton phase diversity using provided Fourier operators
    return +((ift*(mods[i] .* phasor.(ft*(guess .* phis[i]))) .* conj.(phis[i]) for i=1:M)...) ./ M
end

function pdgs(imgs::NTuple{M,LF{Intensity,<:Number,N}},divPhases::NTuple{M,LF{<:Phase,<:Number,N}},nit::Integer,beamGuess::LF{ComplexAmp,<:Complex,N}) where {M,N}
    return pdgs(Tuple(sqrt(i) for i in imgs),divPhases,nit,beamGuess)
end


#---------------------- One shot diversity estimates ------------------------------

function oneShot(img::LF{Intensity,<:Real,N},alpha::Real,beta::NTuple{N,Real},sc::Real) where N
    dL = dualShiftLattice(img.L,img.flambda)
    xc = beta .* (img.flambda / sc)
    divPhase = LF{RealPhase}([alpha/(2*sc^2) * sum(dL[j][I[j]]^2 for j=1:N) + sum(beta[j]*dL[j][I[j]] for j=1:N)/sc for I in CartesianIndices(size(b))])
    dualDivPhase = LF{RealPhase}([ sc^2/(2*alpha*flambda^2) * sum((L[j][I[j]] - xc[j])^2 for j=1:N) for I in CartesianIndices(size(img))], img.L, img.flambda)
    x = sqrt(img) * dualDivPhase
    return LF{ComplexAmp}(x.data |> ifftshift |> ifft |> fftshift,dL,x.flambda) * conj(divPhase)
end