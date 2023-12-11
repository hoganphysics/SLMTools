module VisualizationHelpers
using Images: Gray
using ..LatticeTools
export look

cycle1(x::Complex) = (angle(x) + pi) / (2pi)

"""
    look(x::LF...), look(f::LF{RealPhase}), look(f::LF{ComplexPhase}), look(f::LF{Modulus}), look(f::LF{ComplexAmp}), look(f::LF{Intensity})

    Generate a heatmap of the given attribute of a LatticeField, or a combined heatmap of multiple attributes.

    # Arguments
    - `f::LF{attribute}

    # Returns
    - A heatmap(s) of the input field.

    This function generates a heatmap of attribute of a LatticeField, or a combined heatmap of multiple attributes. The following attributes are supported:
    - `RealPhase`: The real part of the phase of the field.
    - `ComplexPhase`: The phase of the field.
    - `Modulus`: The modulus of the field.
    - `ComplexAmp`: The complex amplitude of the field.
    - `Intensity`: The intensity of the field.
    - all of them together
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