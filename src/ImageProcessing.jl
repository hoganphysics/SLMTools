module ImageProcessing
using FileIO: load
using Images: Gray, RGB

export getImagesAndFilenames, imageToFloatArray, itfa, getCamGrid


"""
    getImagesAndFilenames(dir::String, extension::String) -> (Array{Image}, Array{String})

    Load images from a specified directory `dir` with a given file `extension`.
    # Arguments
    - `dir::String`: The directory from which to load images.
    - `extension::String`: The extension of the image files to be loaded (e.g., ".png", ".jpg").
    # Returns
    - `Array{Image}`: An array of images loaded from the specified directory.
    - `Array{String}`: An array of filenames corresponding to the loaded images.
    """
function getImagesAndFilenames(dir::String, extension::String)
    L = length(extension)                                   # Extension length
    f = s -> length(s) > L && s[end-L+1:end] == extension    # Filter function
    fileNames = filter(f, readdir(dir))
    images = load.(dir .* fileNames)
    return images, fileNames
end

"""
    imageToFloatArray(img::Matrix) -> Matrix{Float64}

    Converts an image to a matrix of Float64, representing grayscale values in the range [0, 1].
    # Arguments
    - `img::Matrix`: The input image matrix, typically in color.
    # Returns
    - `Matrix{Float64}`: The resulting grayscale image as a matrix of floats.
    """
imageToFloatArray(img::Matrix) = convert.(Float64, Gray.(img))

# Alias for `imageToFloatArray`
itfa = imageToFloatArray;


"""
    getCamGrid(linImgs::Vector{Array{T,2}}, idxs, dxcam; roi=nothing)
    Processes a collection of linear phase RGB images to determine the camera grid orientation and origin.

    # Arguments
    - `linImgs::Vector{Array{T,2}}`: A vector of 2D RGB images.
    - `idxs`: Indices or parameters associated with each image in `linImgs`.
    - `dxcam`: The physical unit conversion factor for camera pixels.
    - `roi` (optional): Region of interest specified as a tuple of ranges. If not provided, the entire image is used.

    # Returns
    - A tuple representing the camera grid in physical units, adjusted by the origin.
    - `θ`: The angle of the grid with respect to some reference.
"""
function getCamGrid(linImgs::Vector{Array{T,2}}, idxs, dxcam; roi=nothing) where {T<:Number}
    isnothing(roi) && (roi = ((1:s for s in size(linImgs[1]))...,))
    imgs = [i[roi...] for i in linImgs]
    θ, origin, s = sweepCoords(imgs, idxs; normalize=true)
    # pixels = ((1:length(r) for r in roi)...,)
    # return (((pixels[j] .- origin[j]) .* dxcam for j = 1:length(pixels))...,), θ

    # more julian way of doing the above, which the linter was cranky about 
    # Assuming roi is a tuple of ranges, we obtain the axes directly
    pixels = map(axes(roi)) do r
        1:length(r)  # This is safe because `length(r)` will correspond to the range's size.
    end

    # We can use `eachindex` to iterate over each index of `pixels`
    # Note that `pixels` is a tuple of ranges, so `eachindex` is equivalent to `1:length(pixels)`
    return tuple((map(j -> (pixels[j] .- origin[j]) .* dxcam, eachindex(pixels))...)), θ

end
getCamGrid(linImgs::Vector{Array{T,2}}, idxs, dxcam; roi=nothing) where {T<:RGB} = getCamGrid(itfa.(linImgs), idxs, dxcam; roi=roi)






end # module ImageProcessing