# LatticeFields.jl
# Tools for handling SLM fields defined on lattices.
# Requires: Interpolations

const Lattice{N} = NTuple{N,AbstractRange} where N 
elq(x::Lattice{N},y::Lattice{N}) where N = (all(isapprox.(x,y)) || throw(DomainError((x,y), "Unequal lattices.")); return nothing;)

#-------------- Field value types ------------------------------
abstract type FieldVal end
abstract type Generic <: FieldVal end

# Phases
abstract type Phase <: FieldVal end
abstract type RealPhase <: Phase end
const UPhase = RealPhase
const UnwrappedPhase = RealPhase
abstract type ComplexPhase <: Phase end
const S1Phase = ComplexPhase

# Intensity
abstract type Intensity <: FieldVal end

# Amplitudes 
abstract type Amplitude <: FieldVal end
# Non-negative real amplitude
abstract type Modulus <: Amplitude end
const RealAmplitude = Modulus
const RealAmp = Modulus
# Complex amplitude
abstract type ComplexAmplitude <: Amplitude end
const ComplexAmp = ComplexAmplitude


#------------- LatticeFields ---------------------------------------------------

struct LatticeField{S<:FieldVal,T,N} #<: AbstractArray{T,N}
    data::AbstractArray{T,N}
    L::Lattice{N}
    flambda::Real
    function LatticeField{S,T,N}(data::AbstractArray{T,N},L::Lattice{N},flambda::Real=1.0) where {S<:FieldVal,T,N}
        if size(data) !== length.(L)
            throw(DimensionMismatch("Field data size does not match lattice size."))
        else
            return new(data,L,flambda)
        end
    end
end
LatticeField{S}(array::AbstractArray{T,N},L::Lattice{N},flambda::Real=1.0) where {S<:FieldVal,T,N} = LatticeField{S,T,N}(array,L,flambda)
const LF = LatticeField

# Minimal methods to make these types function as arrays.  
# All arithmetic works on a LatticeField, though the output will be some sort of array, and not a LatticeField. 
Base.size(f::LatticeField) = size(f.data)
Base.getindex(f::LatticeField, i::Int) = getindex(f.data,i)
Base.setindex!(f::LatticeField, v, i::Int) = setindex!(f.data, v, i)
Base.IndexStyle(::Type{<:LatticeField})= IndexLinear()
Base.axes(f::LatticeField,d) = axes(f.data,d)
Base.show(io::IO, f::T) where T<:LF = (println(T); println("Lattice: ",f.L); display(f.data))

#= 
# The following functions allow for creation of sublattice fields by subscripting. But there are some 
# downsides to allowing this; for example, you can get errors by calling functions like sqrt on a lattice field. 
# I've decided for now to leave this unimplmented, and have subfield functions that accomplish this instead. 
Update: I've re-implemented these functions after making LF not a subtype of AbstractArray.
=#
function Base.getindex(f::LatticeField{S,T,N}, x::I...) where {S,T,N,I<:Integer}
    # Subscripting f by integers returns the field value at that location. 
    return getindex(f.data,x...)
end
function Base.getindex(f::LatticeField{S,T,N}, x::AbstractRange{I}...) where {S,T,N,I<:Integer}
    # Subscripting f by ranges returns a subfield on the sublattice. 
    return LatticeField{S,T,N}(getindex(f.data,x...),sublattice(f.L,x...),f.flambda)
end
function Base.getindex(f::LatticeField{S,T,N}, x::Union{I,AbstractRange{I}}...) where {S,T,N,I<:Integer}
    # Subscripting f by integers and ranges returns a subfield on the slice fixing dimensions with integer indices. 
    # Note that this is not a type safe method, so it should not be used in performance demanding situations. 
    z = filter(i->!(x[i] isa Integer), 1:length(x))
    return LatticeField{S,T,length(z)}(getindex(f.data,x...),sublattice(f.L[z],x[z]...),f.flambda)
end 

function subfield(f::LatticeField{S,T,N}, x::I...) where {S,T,N,I<:Integer}
    # Subscripting f by integers returns the field value at that location. 
    return getindex(f.data,x...)
end
function subfield(f::LatticeField{S,T,N}, x::AbstractRange{I}...) where {S,T,N,I<:Integer}
    # Subscripting f by ranges returns a subfield on the sublattice. 
    return LatticeField{S,T,N}(getindex(f.data,x...),sublattice(f.L,x...),f.flambda)
end
function subfield(f::LatticeField{S,T,N}, x::Union{I,AbstractRange{I}}...) where {S,T,N,I<:Integer}
    # Subscripting f by integers and ranges returns a subfield on the slice fixing dimensions with integer indices. 
    # Note that this is not a type safe method, so it should not be used in performance demanding situations. 
    z = filter(i->!(x[i] isa Integer), 1:length(x))
    return LatticeField{S,T,length(z)}(getindex(f.data,x...),sublattice(f.L[z],x[z]...),f.flambda)
end

# Equal lattice query
elq(x::LF,y::LF) = (all(isapprox.(x.L,y.L)) || throw(DomainError((x.L,y.L), "Unequal lattices.")); return nothing;)

#-------------------- Lattice math and conversions -------------------------------

# Default behavior is to throw an undefined method error.
Base.:(*)(args::LF...) = error("Behavior undefined for this combination of inputs: ",*((string(i)*", " for i in typeof.(args))...))
Base.:(+)(args::LF...) = error("Behavior undefined for this combination of inputs: ",*((string(i)*", " for i in typeof.(args))...))
Base.:(-)(args::LF...) = error("Behavior undefined for this combination of inputs.",*((string(i)*", " for i in typeof.(args))...))
Base.:(/)(args::LF...) = error("Behavior undefined for this combination of inputs.",*((string(i)*", " for i in typeof.(args))...))
Base.sqrt(args::LF...) = error("Behavior undefined for this combination of inputs.",*((string(i)*", " for i in typeof.(args))...))
square(args::LF...) = error("Behavior undefined for this combination of inputs.",*((string(i)*", " for i in typeof.(args))...))
Base.abs(args::LF...) = error("Behavior undefined for this combination of inputs.",*((string(i)*", " for i in typeof.(args))...))

# Intensity to Modulus
Base.sqrt(f::LF{Intensity}) = LF{Modulus}(sqrt.(f.data),f.L,f.flambda)

# Amplitude times phase
Base.:(*)(x::LF{<:Amplitude},y::LF{ComplexPhase}) = (elq(x,y); LF{ComplexAmp}(x.data .* y.data, x.L,x.flambda))
Base.:(*)(x::LF{<:Amplitude},y::LF{RealPhase}) = (elq(x,y); LF{ComplexAmp}(x.data .* exp.(2pi*im*y.data), x.L,x.flambda))
Base.:(*)(y::LF{ComplexPhase},x::LF{<:Amplitude}) = (elq(x,y); LF{ComplexAmp}(x.data .* y.data, x.L,x.flambda))
Base.:(*)(y::LF{RealPhase},x::LF{<:Amplitude}) = (elq(x,y); LF{ComplexAmp}(x.data .* exp.(2pi*im*y.data), x.L,x.flambda))

# Combining phases
Base.:(*)(x::LF{ComplexPhase},y::LF{ComplexPhase}) = (elq(x,y); LF{ComplexPhase}(x.data .* y.data, x.L,x.flambda))
Base.:(*)(x::LF{RealPhase},y::LF{ComplexPhase}) = (elq(x,y); LF{ComplexPhase}(exp.(2pi*im*x.data) .* y.data, x.L,x.flambda))
Base.:(*)(x::LF{ComplexPhase},y::LF{RealPhase}) = (elq(x,y); LF{ComplexPhase}(x.data .* exp.(2pi*im*y.data), x.L,x.flambda))

Base.:(+)(x::LF{RealPhase},y::LF{RealPhase}) = (elq(x,y); LF{RealPhase}(x.data .+ y.data, x.L,x.flambda))

# ComplexAmplitude to Modulus
Base.abs(x::LF{ComplexAmp}) = LF{Modulus}(abs.(x.data),x.L,x.flambda)

# Modulus or ComplexModulus to Intensity
square(x::LF{<:Amplitude}) = LF{Intensity}(abs.(x.data).^2,x.L,x.flambda)

# RealPhase to ComplexPhase
wrap(x::LF{RealPhase}) = LF{ComplexPhase}(exp.(2pi*im*x.data),x.L,x.flambda)
wrap(x::LF{ComplexPhase}) = copy(x)

# Phase conjugation
Base.conj(x::LF{RealPhase}) = LF{RealPhase}(-x.data,x.L,x.flambda)
Base.conj(x::LF{ComplexPhase}) = LF{ComplexPhase}(conj.(x.data),x.L,x.flambda)

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

#------------------- Padding -------------------------------------

padout(L::AbstractRange,p::Tuple{Int,Int}) = (0:sum(p)+length(L)-1) .* step(L) .+ (L[1]-step(L)*p[1])
padout(L::AbstractRange,p::Int) = padout(L,(p,p))
padout(L::Lattice{N},p::NTuple{M,Int}) where {N,M} = M==2N ? ( (padout(L[i],(p[2i-1],p[2i])) for i=1:N)..., ) : error("Bad tuple length.")
padout(L::Lattice{N},p::NTuple{N,Int}) where N =  ( (padout(L[i],p[i]) for i=1:N)..., )
padout(L::Lattice{N},p::Int) where N = ( (padout(L[i],p) for i=1:N)..., )

function padout(A::Array{T,N},p::NTuple{M,Int},filler=zero(T)) where {T,N,M}
    M == 2N || error("Bad tuple length.")
    output = fill(filler,size(A) .+ ((p[2i-1]+p[2i] for i=1:N)...,))
    output[((p[2i-1]+1:size(output,i)-p[2i]) for i=1:N)...] .= A
    return output
end
padout(A::Array{T,N},p::NTuple{N,Int},filler=zero(T)) where {T,N} = padout(A,((p[ceil(Int,i/2)] for i=1:2N)...,),filler)
padout(A::Array{T,N},p::Int,filler=zero(T)) where {T,N} = padout(A,((p for i=1:N)...,),filler)

padout(f::LF{S,T,N}, p::NTuple{N,Int},filler=zero(T)) where {S<:FieldVal,T,N} = LF{S,T,N}(padout(f.data,p,filler),padout(f.L,p),f.flambda)


#------------------------ Downsampling and upsampling lattices -------------------------

function sublattice(L::Lattice{N},box::CartesianIndices{N}) where N
    return ((L[j][box[1][j]:box[end][j]] for j=1:N)...,)
end
function sublattice(L::Lattice{N},x::Union{I,AbstractRange{I}}...) where {N,I<:Integer}
    length(L) != length(x) && throw(DimensionMismatch("Wrong number of indices while attempting to make sublattice."))
    y = (( (i isa Integer) ? (i:i) : i for i in x)...,)
    return ((L[j][y[j]] for j=1:N)...,)
end

function upsample(r::AbstractRange,n::Int)
    # Upsamples a range r, replacing each original point by n points, centered on the original. 
    return ((1:length(r)*n) .- (1+n)/2) .* (step(r)/n) .+ r[1]
end
function upsample(L::Lattice{N},n::Int) where N
    # Upsamples a lattice, replacing each original pixel by n^N pixels, centered on the original. 
    return ((upsample(l,n) for l in L)...,)
end
function upsample(L::Lattice{N},ns::NTuple{N,Int}) where N
    # Upsamples a lattice, replacing each original pixel by n^N pixels, centered on the original. 
    return ((upsample(L[i],ns[i]) for i=1:N)...,)
end

function downsample(r::AbstractRange,n::Int)
    # Downsamples a range r by merging n points into one, centered at the mean of the original points.
    # n must divide the length of r. 
    if length(r) % n != 0
        throw(DomainError((n,length(r)), "downsample: Downsample factor n does not divide the length of the range r."))
    end
    return (0:length(r)÷n - 1) .* (step(r)*n) .+ (r[n]+r[1])/2
end
function downsample(L::Lattice{N},n::Int) where N
    # Downsamples a lattice L by merging n^N points into one, centered at the mean of the original points.
    # n must divide the length of each dimension of L. 
    return ((downsample(l,n) for l in L)...,)
end
function downsample(L::Lattice{N},ns::NTuple{N,Int}) where N
    # Downsamples a lattice L by merging ns[1]×ns[2]×⋯×ns[N] points into one, centered at the mean of the original points.
    # ns[i] must divide the length of L[i] for i=1:N. 
    return ((downsample(L[i],ns[i]) for i=1:N)...,)
end

#------------------------ Downsampling and upsampling arrays -------------------------


function upsample(x::AbstractArray{T,N},Ld::Lattice{N},Lu::Lattice{N}; interpolation=cubic_spline_interpolation, bc=Periodic()) where {T<:Number,N}
    # Upsamples an array from coarse lattice Ld to fine lattice Lu.  Actually, there doesn't need to be any relationship
        # between the lattices at all--this function just interpolates x from one lattice to the other. 
    if size(x) != length.(Ld)
        throw(DomainError((size(x),length.(Ld)), "upsample: Size of array x does not match size of lattice Ld."))
    end
    return interpolation(Ld,x,extrapolation_bc=bc)[Lu...]
end
function upsample(x::AbstractArray{T,N},n::Int; interpolation=cubic_spline_interpolation, bc=Periodic()) where {T<:Number,N}
    # Upsamples an array from to a grid that is n times finer, i.e. a grid upsampled by n.
    Ld = ((1:s for s in size(x))...,)
    Lu = upsample(Ld,n)
    return upsample(x,Ld,Lu; interpolation=interpolation, bc=bc)
end
function upsample(x::AbstractArray{T,N},ns::NTuple{N,Int}; interpolation=cubic_spline_interpolation, bc=Periodic()) where {T<:Number,N}
    # Upsamples an array from to a grid that is n times finer, i.e. a grid upsampled by n.
    Ld = ((1:s for s in size(x))...,)
    Lu = upsample(Ld,ns)
    return upsample(x,Ld,Lu; interpolation=interpolation, bc=bc)
end

# This lattice downsampler is basically the same as the upsampler with the roles of Lu and Ld reversed. 
function downsample(x::AbstractArray{T,N},Lu::Lattice{N},Ld::Lattice{N}; interpolation=cubic_spline_interpolation, bc=Periodic()) where {T<:Number,N}
    # Downsamples an array from fine lattice Ld to coarse lattice Lu.  Actually, there doesn't need to be any relationship
        # between the lattices at all--this function just interpolates x from one lattice to the other. 
    if size(x) != length.(Lu)
        throw(DomainError((size(x),length.(Lu)), "upsample: Size of array x does not match size of lattice Lu."))
    end
    return interpolation(Lu,x,extrapolation_bc=bc)[Ld...]
end
function downsample(x::AbstractArray{T,N},n::Int; interpolation=cubic_spline_interpolation, bc=Periodic()) where {T<:Number,N}
    # Downsamples an array from to a grid that is n times coarser, i.e. a grid downsampled by n.
    Lu = ((1:s for s in size(x))...,)
    Ld = downsample(Lu,n)
    return downsample(x,Lu,Ld; interpolation=interpolation, bc=bc)
end
function downsample(x::AbstractArray{T,N},ns::NTuple{N,Int}; interpolation=cubic_spline_interpolation, bc=Periodic()) where {T<:Number,N}
    # Downsamples an array from to a grid that is n times coarser, i.e. a grid downsampled by n.
    Lu = ((1:s for s in size(x))...,)
    Ld = downsample(Lu,ns)
    return downsample(x,Lu,Ld; interpolation=interpolation, bc=bc)
end

# This is also a downsampler, but it averages over superpixels, rather than interpolating. 
function coarsen(x::AbstractArray{T,N},ns::NTuple{N,Int};reducer=(x::AbstractArray -> sum(x)/length(x[:]))) where {T<:Number,N}
    # Downsamples an array to a grid that is ns[i] times coarser in dimension i, i.e. a grid downsampled by ns, 
        # reducing superpixels via the supplied function reducer. 
    sx = size(x)
    if sx .% ns != size(x).*0
        throw(DomainError((sx,ns), "coarsen: Downsample factors ns do not divide the size of array x."))
    end
    box = CartesianIndices(ns) .- CartesianIndex((1 for i=1:N)...)
    CI = CartesianIndices( ((1:ns[i]:size(x)[i] for i=1:N)...,) )
    return [reducer(x[I .+ box]) for I in CI]
end
function coarsen(x::AbstractArray{T,N},n::Int;reducer=(x::AbstractArray -> sum(x)/length(x[:]))) where {T<:Number,N}
    # Downsamples an array to a grid that is n times coarser, i.e. a grid downsampled by n, reducing superpixels via the supplied function reducer. 
    return coarsen(x,((n for i=1:N)...,); reducer=reducer)
end

#------------------------ Downsampling and upsampling fields -------------------------

function upsample(f::LatticeField{S,T,N},Lu::Lattice{N}; interpolation=cubic_spline_interpolation, bc=Periodic()) where {S<:FieldVal,T<:Number,N}
    # Upsamples field f to lattice Lu.  Actually, there doesn't need to be any relationship
        # between the lattices at all--this function just interpolates x from one lattice to the other. 
    return LatticeField{S}(upsample(f.data,f.L,Lu; interpolation=interpolation, bc=bc),Lu,f.flambda)
end
function upsample(f::LatticeField{S,T,N},n::Int; interpolation=cubic_spline_interpolation, bc=Periodic()) where {S<:FieldVal,T<:Number,N}
    # Upsamples field f to a grid that is n times finer, i.e. a grid upsampled by n.
    Lu = upsample(f.L,n)
    return upsample(f,Lu; interpolation=interpolation, bc=bc)
end
function upsample(f::LatticeField{S,T,N},ns::NTuple{N,Int}; interpolation=cubic_spline_interpolation, bc=Periodic()) where {S<:FieldVal,T<:Number,N}
    # Upsamples an array from to a grid that is ns[i] times finer in dimension i, i.e. a grid upsampled by n.
    Lu = upsample(f.L,ns)
    return upsample(f,Lu; interpolation=interpolation, bc=bc)
end

function downsample(f::LatticeField{S,T,N},Ld::Lattice{N}; interpolation=cubic_spline_interpolation, bc=Periodic()) where {S<:FieldVal,T<:Number,N}
    # Downsamples field f to lattice Lu.  Actually, there doesn't need to be any relationship
        # between the lattices at all--this function just interpolates x from one lattice to the other. 
    return LatticeField{S}(downsample(f.data,f.L,Ld;interpolation=interpolation,bc=bc),Ld,f.flambda)
end
function downsample(f::LatticeField{S,T,N},n::Int; interpolation=cubic_spline_interpolation, bc=Periodic()) where {S<:FieldVal,T<:Number,N}
    # Downsamples field f to a grid that is n times coarser, i.e. a grid downsampled by n.
    Ld = downsample(f.L,n)
    return LatticeField{S}(downsample(f.data,f.L,Ld;interpolation=interpolation,bc=bc),Ld,f.flambda)
end
function downsample(f::LatticeField{S,T,N},ns::NTuple{N,Int}; interpolation=cubic_spline_interpolation, bc=Periodic()) where {S<:FieldVal,T<:Number,N}
    # Downsamples an array to a grid that is ns[i] times coarser in dimension i, i.e. a grid downsampled by n.
    Ld = downsample(f.L,ns)
    return LatticeField{S}(downsample(f.data,f.L,Ld;interpolation=interpolation,bc=bc),Ld,f.flambda)
end

function coarsen(f::LatticeField{S,T,N},ns::NTuple{N,Int};reducer=(x::AbstractArray -> sum(x)/length(x[:]))) where {S<:FieldVal,T<:Number,N}
    # Downsamples field f to a grid that is ns[i] times coarser in dimension i, i.e. a grid downsampled by ns, 
        # reducing superpixels via the supplied function reducer. 
    Ld = downsample(f.L,ns)
    return LatticeField{S}(coarsen(f.data,ns;reducer=reducer),Ld,f.flambda)
end
function coarsen(f::LatticeField{S,T,N},n::Int;reducer=(x::AbstractArray -> sum(x)/length(x[:]))) where {S<:FieldVal,T<:Number,N}
    # Downsamples an array from to a grid that is n times coarser, i.e. a grid downsampled by n, reducing superpixels via the supplied function reducer. 
    Ld = downsample(f.L,n)
    return LatticeField{S}(coarsen(f.data,n;reducer=reducer),Ld,f.flambda)
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

