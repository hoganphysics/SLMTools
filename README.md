# SLMTools

#### Created by Hogan Lab, Stanford University.

This package is designed to facilitate working with a spatial light modulator (SLM) in the context of laser beam shaping.  
It provides various algorithms for calibrating input beams and generating SLM phases. 

## Example

```julia
using SLMTools

# Generate an input grid and corresponding output grid
N = 128
L0 = natlat((N,N))
dL0 = dualShiftLattice(L0)

# Generate an input beam and target output beam
inputBeam = lfGaussian(Intensity, L0, 1.0)
targetBeam = lfRing(Intensity, dL0, 2.5, 0.5)

# Use optimal transport to find an SLM phase to make an approximate output beam
phiOT = otPhase(inputBeam,targetBeam,0.001)

# Refine the above phase using the Gerchberg-Saxton algorithm.
phiGS = gs(inputBeam,targetBeam,100,phiOT)

# View the resulting output beams
outputOT = square(sft(sqrt(inputBeam) * phiOT))
outputGS = square(sft(sqrt(inputBeam) * phiGS))
look(targetBeam,outputOT,outputGS)
```

## Installation
Clone this repository into a local directory `dir/to/repo`.  From the Julia REPL, run `]add dir/to/repo`. 

## Package overview
The basic object in this package is the `LatticeField`, or `LF` for short.  An LF is a representation of a field sampled on a rectangular grid.  It includes the field values at each point of the grid (represented as an array), the grid itself, and a scale factor which is important for modelling SLMs, namely the product of the light wavelength and a lens focal length.  The basic definition is as follows: 
```
struct LatticeField{S<:FieldVal,T,N}
    data::AbstractArray{T,N}
    L::Lattice{N}
    flambda::Real
end
```
The full definition of `LF` is a bit more involved to enforce some constraints, but this is the main idea.  The `Lattice{N}` type used above is just another name for a tuple of `AbstractRange` objects, and represents the grid.  Here's the full definition: 
```
const Lattice{N} = NTuple{N,AbstractRange} where {N}
```
The `S<:FieldVal` type parameter in the definition of `LF` is a way to keep track of what the particular `LF` represents, e.g. an `Intensity`, `Phase`, `Modulus`, etc.  This provides very useful error checking on your code.  The currently supported `FieldVal` subtypes are: 
* `Generic`
* `Phase`
* `RealPhase <: Phase` (aliases: `UPhase`, `UnwrappedPhase`)
* `ComplexPhase <: Phase` (alias: `S1Phase`)
* `Intensity`
* `Amplitude`
* `Modulus <: Amplitude` (alias: `RealAmplitude`, `RealAmp`)
* `ComplexAmplitude <: Amplitude` (alias: `ComplexAmp`)
* 
These types don't do anything on their own--they are only used for keeping track of what each `LF` represents and defining behavior of function methods (i.e. multiple dispatch). 

The basic workflow of SLMTools is this:
* Define two `Lattice` objects representing the pixel grid for your SLM and camera. 
* Define some `LF` objects representing phases, intensities, etc. in the camera and SLM plane.
* Use functions like `otPhase` or `gs` to determine what SLM phases will shape a given input beam into a target output beam.
* Use functions like `pdgs` to infer the input beam on your SLM from intensities measured on your camera.

### Useful functions
The following are some of the major useful functions provided by this package. 
* `otPhase`: Uses optimal transport to generate an SLM phase for transforming a given input beam into a target output beam.
* `gs`: Uses the Gerchberg-Saxton algorithm to generate an SLM phase for transforming a given input beam into a target output beam.
* `mraf`: Uses the Mixed Region Amplitude Freedom (MRAF) algorithm to generate an SLM phase for transforming a given input beam into a target output beam.
* `pdotBeamEstimate`: Uses optimal transport to estimate the input beam incident upon an SLM, given a set of diversity phase images and phase coefficients.
* `pdgs`: Uses an iterative Fourier transform algorithm to estimate the input beam incident upon an SLM, given a set of diversity phase images and phase coefficients.

Most functions have doc strings detailing usage, and you should refer to those for more info. 
