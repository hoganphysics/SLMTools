module HoganParameters
using LatticeTools: dualShiftLattice, padout

"""
    dxcam
    camera pixel size in microns
"""
dxcam = 5.86  

"""
    dxslm
    SLM pixel size in microns
"""
dxslm = 17

"""
    nslm
    SLM pixel count, assuming square SLM
"""
nslm = 1024

"""
    flambda
    wavelength of light in ... tens of picometers?
"""
flambda = 1.064 * 10^5

"""
    dXslm
    the SLM pixel size in something?
"""
dXslm = flambda / (nslm * dxslm)

"""
    Lslm
    SLM lattice in microns, centered at (0,0)
"""
Lslm = ((1:1024) .- 513, (1:1024) .- 513) .* 17

"""
    dL
    dual lattice of Lslm
"""
dL = dualShiftLattice(Lslm, flambda)

"""
    LslmE
    padded Lslm
"""
LslmE = padout(Lslm, 100)

"""
    dLE
    dual lattice of LslmE
"""
dLE = dualShiftLattice(LslmE, flambda)


end # module HoganParameters