module ImageProcessing
using FileIO: load
using Images: Gray

export getImagesAndFilenames, imageToFloatArray, itfa


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







end # module ImageProcessing