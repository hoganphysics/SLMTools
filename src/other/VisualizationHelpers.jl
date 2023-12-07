module VisualizationHelpers
using Plots: heatmap, plot
"""
    hm(a)

    Create a heatmap of the input array with no colorbar.

    # Arguments
    - `a`: An array for which the heatmap is to be generated.

    # Returns
    - A heatmap of the input array `a`.

    # Overloaded Method
    - `hm(a...)`: Generates a combined plot of heatmaps for each input array in `a`.

    This function and its overloaded method are convenience wrappers for creating heatmaps without color bars, useful for visualizing matrices or arrays of data.
"""
hm(a) = heatmap(a, colorbar=false)
hm(a...) = plot((hm(i) for i in a)...)

"""
    plotBeamMap(beams::Vector{M}; st=1, roi=nothing) where {M<:Matrix}

    Generate a series of heatmaps for the intensity and phase of each beam in a vector of beams.

    # Arguments
    - `beams::Vector{M}`: A vector where each element is a matrix representing a complex beam.
    - `st`: Step size for the region of interest (ROI) in the plot, defaults to 1.
    - `roi`: An optional tuple specifying the ROI over which to plot the beams.

    # Returns
    - A combined plot of heatmaps showing the normalized absolute square and angle of each beam within the specified ROI.

    This function iterates over the vector of beams, creating two heatmaps for each: one for the normalized intensity (absolute square) and one for the phase (angle) of the beam, then vertically concatenates them for a comprehensive visualization.
"""
function plotBeamMap(beams::Vector{M}; st=1, roi=nothing) where {M<:Matrix}
    if isnothing(roi)
        roi = (1:st:s for s in size(beams[1]))
    end
    plot(vcat(([heatmap(nabs(beam[roi...]) .^ 2, colorbar=false), heatmap(angle.(beam[roi...]), colorbar=false)] for beam in beams)...)...)
end


"""
    plotBeamMap(beam; st=1, roi=nothing)

    Generate heatmaps for the intensity and phase of a single beam.

    # Arguments
    - `beam`: A matrix representing a complex beam.
    - `st`: Step size for the region of interest (ROI) in the plot, defaults to 1.
    - `roi`: An optional tuple specifying the ROI over which to plot the beam.

    # Returns
    - A plot with two heatmaps: one showing the normalized absolute square and the other showing the angle of the beam within the specified ROI.

    This function creates a heatmap for the normalized intensity (absolute square) and another for the phase (angle) of the given beam and displays them side by side.
"""
function plotBeamMap(beam; st=1, roi=nothing)
    if isnothing(roi)
        roi = (1:st:s for s in size(beam))
    end
    plot(heatmap(nabs(beam[roi...]) .^ 2, colorbar=false), heatmap(angle.(beam[roi...]), colorbar=false))
end

end # module VisualizationHelpers