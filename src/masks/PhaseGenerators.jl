module PhaseGenerators
using FileIO: save
using Images: Gray, n0f8

export makeRampedParabola
"""
    makeRampedParabola(imSize::Tuple{Integer,Integer}, pixelScale::Real, center::Tuple{Real,Real}, pAmp::Real, rAmps::Tuple{Real,Real}; outputname=nothing)
Create a grayscale image of a ramped parabolic phase map.

# Arguments
    - `imSize::Tuple{Integer,Integer}`: The size of the image (number of pixels in the y and x dimensions).
    - `pixelScale::Real`: The size of a superpixel in pixels. If pAmp=1, this is the number of pixels from the center to first phase rollover
    - `center::Tuple{Real,Real}`: The center of the parabola in pixel coordinates.
    - `pAmp::Real`: The amplitude of the parabola. An unwrapped phase map will go from 0 to pAmp over pixelScale pixels (ignoring the linear ramp).
    - `rAmps::Tuple{Real,Real}`: The amplitudes of the linear ramp in the y and x dimensions.

    # Keyword Arguments
    - `outputname`: The optional filename to save the generated image. If not provided, the image is returned instead.

    # Returns
    - A grayscale image array representing the ramped parabolic phase map if `outputname` is `nothing`, or `nothing` if the image is saved to a file.

    This function generates a 2D parabolic phase map with a linear ramp, normalizes the phase values to lie within the range [0, 1], and optionally saves the image if `outputname` is provided. The phase map is useful for optical applications where a quadratic phase profile is required, such as lens simulations.   
"""
function makeRampedParabola(imSize::Tuple{Integer,Integer}, pixelScale::Real, center::Tuple{Real,Real}, pAmp::Real, rAmps::Tuple{Real,Real}; outputname=nothing)
    # This should be modified to have a 1/2 in the exponent at some point. 
    img = ([pAmp * sum((Tuple(I) .- center) .^ 2) / pixelScale^2 for I in CartesianIndices(imSize)] .+
          [sum((Tuple(I) .- center) .* rAmps) / pixelScale for I in CartesianIndices(imSize)]) |> (x -> mod.(x, 1))
    if isnothing(outputname)
        return n0f8.(Gray.(img))
    else
        save(outputname, n0f8.(Gray.(img)))
        return nothing
    end
end


end # module PhaseGenerators