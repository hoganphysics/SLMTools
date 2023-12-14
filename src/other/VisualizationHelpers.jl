module VisualizationHelpers
using Images: Gray
using ..LatticeTools
export look

cycle1(x::Complex) = (angle(x) + pi) / (2pi)

"""
    look(x::LF...), look(f::LF{RealPhase}), look(f::LF{ComplexPhase}), look(f::LF{Modulus}), look(f::LF{ComplexAmp}), look(f::LF{Intensity})

    Generate a visualization of a LatticeField, or a combined visualization of multiple attributes.  The flavor of visualization is tailored to the FieldVal type.  

    # Arguments
    - `f::LF{FieldVal}

    # Returns
    - An image of the LatticeField.

    This function generates a heatmap of attribute of a LatticeField, or a combined heatmap of multiple attributes. The following attributes are supported:
    - `RealPhase`: Unwrapped phase is displayed mod 1, with black coloration representing phases just above 0 and white representing phases just below 1.
    - `ComplexPhase`: Wrapped phase is displayed mod 2pi, black representing phases just above 0 and white representing phases just below 2pi.
    - `Modulus`: Modulus is displayed normalized so that the maximum value is fully saturated.
    - `ComplexAmp`: ComplexAmp is displayed as a side-by-side plot of Modulus and Phase, with each individually displayed as above.
    - `Intensity`: Intensity is displayed normalized so that the maximum value is fully saturated.
    - Multiple inputs are displayed as the concatenation of their individual outputs. 
    """
function look(f::LF{RealPhase})
    return Gray.((f.data .- minimum(f.data)) .% 1)
end
function look(f::LF{ComplexPhase})
    return Gray.(cycle1.(f.data))
end
function look(f::LF{Modulus})
    return Gray.(f.data ./ maximum(f.data))
end
function look(f::LF{ComplexAmp})
    return hcat(Gray.(abs.(f.data) ./ maximum(abs.(f.data))), Gray.(cycle1.(f.data)))
end
function look(f::LF{Intensity})
    return Gray.(f.data ./ maximum(f.data))
end
function look(x::LF...)
    return hcat((look(y) for y in x)...)
end

end # module VisualizationHelpers