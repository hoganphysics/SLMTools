

#--------------------- Basic fitting ----------------------------------
"""
    findCOM(img, threshold=0.1)
    Finds the center of mass (COM) of an image based on pixel intensity, after applying a threshold.
    Parameters:
    - img: The input image as a 2D array where each element represents pixel intensity.
    - threshold=0.1: A pixel intensity threshold. Pixels with intensity below this value are not considered in the COM calculation.    #
    Returns:
    - A 2-element array representing the coordinates of the COM of the image.
    Apply the threshold to the image, setting values below the threshold to zero.
"""
function findCOM(img, threshold=0.1)
    clipped = (x -> (x > threshold ? x : zero(x))).(img)
    s = sum(clipped)
    s > 0 ? (clipped /= s) : error("Black image. Can't normalize")
    return sum(I -> clipped[I] * [I[1], I[2]], CartesianIndices(size(clipped)))
end

"""
    lineFit(arr, xs)
    Performs a linear fit to a set of data points.
    Parameters:
    - arr: A 1D array of dependent variable values (typically 'y' values).
    - xs: A 1D array of independent variable values (typically 'x' values) corresponding to 'arr'.    #
    Returns:
    - A tuple containing the slope and intercept of the fitted line.
    Perform the linear fit by solving the normal equations.
    `hcat(xs, ones(length(xs)))` creates a matrix with the independent variable 'xs' and a column of ones for the intercept term.
    The `\` operator performs the matrix division, effectively solving the normal equations for the best fit.
    The result is unpacked into a tuple containing the slope and intercept.
    """
function lineFit(arr, xs)
    return ((hcat(xs, ones(length(xs))) \ arr)...,)
end

"""
    sweepCoords(imgs::Vector{Matrix{T}}, idxs; normalize=false, threshold=0.1, plt=false) where {T<:Number}

    Analyzes a set of images represented as matrices to compute the slope, intercept, 
    and angle of the beam center of mass relative to provided indices or parameters.

    # Arguments
    - `imgs::Vector{Matrix{T}}`: A vector of matrices, where each matrix represents an image.
    - `idxs`: A vector of indices or parameters for analysis.
    - `normalize::Bool` (optional): If true, normalizes each image in `imgs` by its maximum value.
    - `threshold::Float64` (optional): Threshold value used in the `findCOM` function.
    - `plt::Bool` (optional): If true, plots the results of the analysis.

    # Returns
    - `θ`: The angle of the beam center of mass trajectory.
    - `[xi, yi]`: Intercepts of the linear fit on x and y coordinates.
    - `[xs, ys]`: Slopes of the linear fit on x and y coordinates.
"""
function sweepCoords(imgs::Vector{Matrix{T}}, idxs; normalize=false, threshold=0.1, plt=false) where {T<:Number}
    if normalize
        imgs = [i / maximum(i) for i in imgs]
    end
    cs = hcat((findCOM(i, threshold) for i in imgs)...)'# Compute center of mass for each image
    xs, xi, ys, yi = lineFit(cs[:, 1], idxs)..., lineFit(cs[:, 2], idxs)...# Perform linear fit on x and y coordinates separately
    θ = angle(xs + im * ys)
    if plt
        p = scatter(idxs, [cs[:, 1], cs[:, 2]], label=["Center X" "Center Y"])
        plot!(p, idxs, [idxs .* xs .+ xi, idxs .* ys .+ yi], label=["Fit X" "Fit Y"])
        display(p)
    end
    return θ, [xi, yi], [xs, ys]
end

"""
    sweepCoords(imgs::Vector{Matrix{T}}, idxs; kwargs...) where {T<:RGB}
    This version of the function is specifically for images of type `RGB`. It first converts the RGB images to 
    a format suitable for analysis (`itfa=imageToFloatArray`) and then calls the main `sweepCoords` function.
"""
function sweepCoords(imgs::Vector{Matrix{T}}, idxs; kwargs...) where {T<:RGB}
    return sweepCoords(itfa.(imgs), idxs; kwargs...)
end
"""
    sweepCoords(imgs::Vector{Matrix{T}}, idxs; kwargs...) where {T<:RGB}
    Analyzes a set of images with respect to a specified grid to compute the angle, adjusted intercepts, 
    and slopes of the beam center of mass.

    This version extends the functionality to handle a grid of coordinates. It adjusts the intercepts and slopes 
    based on the step sizes of the grid ranges.
"""
function sweepCoords(imgs, idxs, grid::NTuple{2,AbstractRange}; kwargs...)
    θ, i, s = sweepCoords(itfa.(imgs), idxs; kwargs...)
    return θ, (r -> r[1] - step(r)).(grid) .+ step.(grid) .* i, s .* step.(grid)
end

# The following is a useful test of how accurate your linear sweep fit is:
# Comparing expected and measured displacements
#beta1 = 8/1024/dxslm    # In μm, based on the index 1 linear phase having 8 periods over a 1024 pixel SLM. 
#slopePred = beta1 * flambda
#slopeMeas = sqrt(xs^2 + ys^2) * dxcam
#slopePred/slopeMeas - 1

