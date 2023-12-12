# slmIOHelpers.jl
# File input and output functions for working with SLMs and cameras. 
# Requires: Images, FileIO, LatticeFields.jl





#--------------------------- Dualization --------------------------------------------

function dualate(f::LF{S,T,2},L::Lattice{2},center::Vector{<:Real},θ::Real,flambda::Real=1.0; 
        roi=nothing, interpolation=cubic_spline_interpolation, naturalize=false,bc=zero(T)) where {S<:FieldVal,T<:Number}
    # Corrects for offset and tilt of f.L and then interpolates f onto a lattice dual to L. 
    length(center) == length(f.L) || error("Incompatible lengths for center and f.L.")
    if !isnothing(roi)
        f = f[roi...]
    end
    L0 = Tuple(f.L[i] .- center[i] for i=1:2)
    interp = interpolation(L0, f.data, extrapolation_bc=bc)
    dL = dualShiftLattice(L,flambda)
    Lx0, Lxs, Ly0, Lys = dL[1][1], step(dL[1]), dL[2][1], step(dL[2])
    R = [cos(θ) -sin(θ) ; sin(θ) cos(θ)]
    r0, dx, dy = (R * [Lx0, Ly0]), Lxs * R[:,1], Lys * R[:,2]
    if naturalize
        return LF{S}([interp((r0 + dx*I[1] + dy*I[2])...) for I in CartesianIndices(length.(dL))],natlat(length.(dL)),1.0);
    else
        return LF{S}([interp((r0 + dx*I[1] + dy*I[2])...) for I in CartesianIndices(length.(dL))],dL,flambda);
    end
end
dualate(fs::Vector{<:LF{S}},L::Lattice{2},center::Vector{<:Real},theta::Real,flambda::Real=1.0;kwargs...) where {S<:FieldVal} = 
    [dualate(f,L,center,theta,flambda;kwargs...) for f in fs]

#--------------------------- File output --------------------------------------------

function savePhase(x::AbstractArray{ComplexF64,2},name::String)
    save(name,Gray.((angle.(x) .+ pi) ./ (2pi)))
end
function savePhase(x::LF{ComplexPhase},name::String)
    savePhase(x.data,name)
end
function savePhase(x::LF{RealPhase},name::String)
    save(name,Gray.(mod.(x.data,1)))
end

#------------------- Shifted Fourier transforms and natlats -------------------------------


function naturalize(f::LF{S,T,N}) where {S<:FieldVal,T,N}
    # Converts LatticeField to be difined on the natural lattice of the same size. 
    return LF{S}(f.data,natlat(size(f)),1.0)
end

#------------------ Standard phases and intensities --------------------------------------


function lfRampedParabola(T::DataType,imSize::NTuple{N,Integer}, pixelScale::Real, center::NTuple{N,Real}, pAmp::Real, rAmps::NTuple{N,Real};L::Union{Nothing,Lattice{N}}=nothing,flambda::Real=1) where N
    if isnothing(L)
        L = natlat(imSize)
    end
    im = ([ pAmp * sum( (Tuple(I) .- center).^2 ) / (2*pixelScale^2) for I in CartesianIndices(imSize)] .+ 
          [ sum( (Tuple(I) .- center) .* rAmps ) / pixelScale  for I in CartesianIndices(imSize)]) |> (x->mod.(x,1))
    return LF{T}(im,L,flambda)
end


function lfGaussian(T::DataType,sz::NTuple{N,Integer},sigma::Number;L::Union{Nothing,Lattice{N}}=nothing,flambda::Real=1) where N
    if isnothing(L)
        L = natlat(sz)
    end
    return LF{T}([exp(-(sum(L[i][I[i]]^2 for i=1:N)/(2*sigma^2))) for I in CartesianIndices(length.(L))],L,flambda)
end

function lfRing(T::DataType,sz::NTuple{N,Integer},r::Number,w::Number;L::Union{Nothing,Lattice{N}}=nothing,flambda::Real=1) where N
    if isnothing(L)
        L = natlat(sz)
    end
    return LF{T}([exp(-(sqrt(sum(L[i][I[i]]^2 for i=1:N))-r)^2/(2*w^2)) for I in CartesianIndices(length.(L))],L,flambda)
end

