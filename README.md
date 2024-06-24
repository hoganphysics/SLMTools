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
inputBeam = lfGaussian(Intensity, (N,N), 1.0 ; L=L0)
targetBeam = lfRing(Intensity, (N,N), 2.5, 0.5; L=L0)

# Use optimal transport to find an SLM phase to make an approximate output beam
phiOT = otPhase(inputBeam,targetBeam,0.001)

# Refine the above phase using the Gerchberg-Saxton algorithm.
phiGS = gs(inputBeam,targetBeam,100,phiOT)

# View the resulting output beams
outputOT = square(sft(sqrt(inputBeam) * phiOT))
outputGS = square(sft(sqrt(inputBeam) * phiGS))
look(targetBeam,outputOT,outputGS)
```

Please see the documentation pages for further information.

## Installation
Clone this repository into a local directory `dir/to/repo`.  From the Julia REPL, run `]add dir/to/repo`. 
