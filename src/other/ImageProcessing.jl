module ImageProcessing
using FileIO: load
using Images: Gray, RGB, Colorant
using ..LatticeTools

export getImagesAndFilenames, imageToFloatArray, itfa, getCamGrid, sweepCoords, findCOM, lineFit, castImage, loadDir, parseFileName,parseStringToNum, getOrientation

"""
    castImage(T::DataType, img::Matrix{<:Colorant}, L::Lattice{2}, flambda::Real)

    Cast an image of colorants into a specified data type and create a LatticeField representation.

    This function takes an image matrix `img`, which consists of colorants (such as RGB values), and casts it 
    into a specific data type `T`. It then creates and returns a `LatticeField` object of type `T` based on the 
    transformed image, the specified lattice `L`, and a wavelength parameter `flambda`.

    # Arguments
    - `T::DataType`: The data type to which the image should be cast.
    - `img::Matrix{<:Colorant}`: A matrix representing the image, where each element is a colorant (e.g., an RGB value).
    - `L::Lattice{2}`: A 2-dimensional lattice that defines the spatial properties of the resulting LatticeField.
    - `flambda::Real`: A real number representing the wavelength parameter for the LatticeField.

    # Returns
    - `LF{T}`: A LatticeField of type `T`, representing the transformed image.

    # Example
    ```julia
    julia> img = load("path/to/image.png") # Assume this returns a Matrix{RGB}
    julia> L = Lattice((10, 10)) # Define a 2D lattice
    julia> flambda = 500.0 # Wavelength parameter
    julia> transformed_img = castImage(Float64, img, L, flambda)
    """
function castImage(T::DataType, img::Matrix{<:Colorant}, L::Lattice{2}, flambda::Real)
    return LF{T}(itfa(img), L, flambda)
end


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
    loadDir(dir::String, extension::String; T::DataType=Intensity, outType::DataType=Float64, 
            L::Union{Lattice{2},Nothing,Real,Tuple{Real,Real}}=nothing, flambda::Real=1.0, 
            cue::Union{String,Char,Nothing}=nothing, look::Union{Symbol,String,Char}=:after)

    Load a directory of images and convert them into LatticeField structures, also parsing their filenames for parameters.

    This function loads all images from a specified directory `dir` with a given file extension `extension`. It converts 
    each image into a LatticeField of type `outType` and optionally extracts parameters from the filenames based on `cue` and `look`.

    # Arguments
    - `dir::String`: The directory from which to load images.
    - `extension::String`: The file extension of images to load (e.g., ".png", ".jpg").
    - `T::DataType`: The data type of the original image data (default: `Intensity`).
    - `outType::DataType`: The output data type of the LatticeFields (default: `Float64`).
    - `L::Union{Lattice{2},Nothing,Real,Tuple{Real,Real}}`: Lattice defining spatial properties. If `Real` or `Tuple`, it's interpreted as pixel spacing.
    - `flambda::Real`: Wavelength parameter for the LatticeFields (default: 1.0).
    - `cue::Union{String,Char,Nothing}`: Optional string or character used as a cue to parse filenames (default: `nothing`).
    - `look::Union{Symbol,String,Char}`: Determines where to look for `cue` in the filename (default: `:after`).

    # Returns
    - A tuple containing:
    - An array of LatticeFields, each representing an image.
    - An array of parameters extracted from the filenames.

    The function ensures all images have the same size and converts them to LatticeFields using `castImage`. 
    Filename parsing is done using `parseFileName`.
    """
function loadDir(dir::String, extension::String;
    T::DataType=Intensity, outType::DataType=Float64, L::Union{Lattice{2},Nothing,Real,Tuple{Real,Real}}=nothing, flambda::Real=1.0,
    cue::Union{String,Char,Nothing}=nothing, look::Union{Symbol,String,Char}=:after)
    # Loads a directory of images, making them into LF{outType} structures, and parses their filenames.  Returns [images], [parameters from filenames].
    images, fileNames = getImagesAndFilenames(dir, extension)
    sz = size(images[1])
    all(size(i) == sz for i in images) || error("Inconsistent image size.")
    if L isa Union{Real,Tuple{Real,Real}}   # Interpret L as the pixel spacing, possibly different in different directions in the case it's a Tuple.
        L = Tuple(1:s for s in sz) .* L
    elseif isnothing(L)
        L = Tuple(1:s for s in sz)
    end
    fields = [castImage(T, i, L, flambda) for i in images]
    if isnothing(cue)
        params = [parseFileName(f; outType=outType) for f in fileNames]
    else
        params = [parseFileName(f, cue, look; outType=outType) for f in fileNames]
    end
    return fields, params
end

"""
    parseStringToNum(s::String; outType=nothing)

    Convert a numeric string to a number, with an option to specify the output type.

    This function takes a string `s` containing numeric characters and converts it into a number. 
    If `s` contains commas, they are treated as decimal points. By default, strings with commas are 
    converted to `Float64`, and those without are converted to `Int`. An optional `outType` can be 
    specified to override this behavior.

    # Arguments
    - `s::String`: A string representing a numeric value. The string should only contain numerals and, optionally, commas.
    - `outType`: An optional data type to which the string should be converted. If not specified, the type is inferred based on the presence of commas.

    # Returns
    - A number of type `Int`, `Float64`, or `outType` (if specified), representing the numeric value of the input string.

    # Notes
    - Commas in the string are replaced with periods to conform with the standard decimal notation.

    The function provides a convenient way to parse numeric strings from various formats, especially useful in data processing where numeric values might be represented as strings with different formatting.
    """
function parseStringToNum(s::String; outType=nothing)
    # Converts a numeric string to a number.  Numeric strings can only contain numerals and commas. 
    # By default, if there is a comma, then the string is converted to a float.  Otherwise it 
    # is converted to an integer. 
    if isnothing(outType)
        if ',' in s
            s = replace(s, ',' => '.')
            return parse(Float64, s)
        else
            return parse(Int, s)
        end
    else
        s = replace(s, ',' => '.')
        return parse(outType, s)
    end
end

"""
    parseFileName(name::String, cue::Union{String,Char}, 
                  look::Union{Symbol,String,Char}=:after; outType=nothing)

    Extract and convert a numeric string from a filename based on a specified cue.

    This function locates a numeric string in `name` based on its relation to a specified `cue` substring or character.
    The position of the numeric string relative to the cue is determined by the `look` parameter.

    # Arguments
    - `name::String`: The filename to parse.
    - `cue::Union{String,Char}`: The substring or character that serves as a reference point in the filename.
    - `look::Union{Symbol,String,Char}`: Specifies whether the numeric string is before or after the cue (default is `:after`).
    - `outType`: An optional data type to which the numeric string should be converted.

    # Returns
    - A number representing the numeric portion of the filename. The type of the number is determined by `outType`, or by the presence of a comma in the numeric string.

    # Notes
    - `look` can be any of `:before`, `:b`, `"before"`, `"b"`, `'b'`, `:after`, `:a`, `"after"`, `"a"`, `'a'`.
    - If `look` starts with 'b', the numeric string is assumed to be before the cue; if with 'a', after the cue.
    - The function assumes that numeric strings can contain only numerals and possibly a comma as a decimal point.
    """
function parseFileName(name::String, cue::Union{String,Char}, look::Union{Symbol,String,Char}=:after; outType=nothing)
    # Finds a numeric string in a filename and turns it into a number.  The numeric string must either precede or follow a 
    # substring or character called "cue".  "look" indicates whether the string follows or precedes the cue.
    # look may be any of :before,:b,"before","b",'b',:after,:a,"after","a",'a'.  All the values starting with a indicate 
    # that the string is after the cue. 
    # The numeric string can only contain numerals and possibly a comma, which stands in for a decimal point.
    # By default, if there's a comma the numeric string is converted to a float, and otherwise it is converted
    # to an integer.  
    look in [:before, :b, "before", "b", 'b', :after, :a, "after", "a", 'a'] || throw(ValueError("Unrecognized look value."))
    look in [:before, :b, "before", "b", 'b'] ? (look = :b) : (look = :a)
    ext = findlast('.', name)
    name = '.' * name[1:ext-1] * '.'    # The periods are to simplify finding the beginning/end of the numeric string. 
    idx = findlast(cue, name)
    if look == :a
        j = findfirst((x -> !(x in [" 0123456789"..., ','])), name[idx[end]+1:end]) + idx[end]
        s = name[idx[end]+1:j-1]
        return parseStringToNum(s, outType=outType)
    else
        j = findlast((x -> !(x in [" 0123456789"..., ','])), name[1:idx[1]-1])
        s = name[j+1:idx[1]-1]
        return parseStringToNum(s, outType=outType)
    end
end

"""
    parseFileName(name::String; outType=nothing)

    Extract and convert the leading numeric string from a filename to a number.

    This function assumes the filename consists of a numeric string followed by a file extension and extracts this numeric portion.

    # Arguments
    - `name::String`: The filename to parse, which should start with a numeric string.
    - `outType`: An optional data type to which the numeric string should be converted.

    # Returns
    - A number representing the leading numeric portion of the filename. The type of the number is determined by `outType`, or by the presence of a comma in the numeric string.

    # Notes
    - The function assumes that numeric strings can contain only numerals and possibly a comma as a decimal point.
    """
function parseFileName(name::String; outType=nothing)
    # Parses a filename to a number.  Assumes that the filename consists of a numeric string followed by a file extension.
    idx = findfirst('.', name)
    return parseStringToNum(name[1:idx-1], outType=outType)
end


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
    # Processes a collection of linear phase images to get the angle and origin of the camera.
    # Outputs the camera grid in physical units and the angle. 
    isnothing(roi) && (roi = ((1:s for s in size(linImgs[1]))...,))
    imgs = [i[roi...] for i in linImgs]
    θ, origin, s = sweepCoords(imgs, idxs; normalize=true)
    pixels = ((1:length(r) for r in roi)...,)
    return (((pixels[j] .- origin[j]) .* dxcam for j = 1:length(pixels))...,), θ
end
getCamGrid(linImgs::Vector{Array{T,2}}, idxs, dxcam; roi=nothing) where {T<:RGB} = getCamGrid(itfa.(linImgs), idxs, dxcam; roi=roi)

"""
    getOrientation(linImgs::Vector{LF{Intensity,T,2}}, idxs::Vector{<:Number}; 
                   roi=nothing, threshold::Number=0.1) where T<:Number

    Determine the orientation of a series of 2D intensity images.

    This function calculates the orientation of a sequence of 2D LatticeFields based on their centroids. 
    It optionally considers a region of interest (ROI) and clips intensities below a specified threshold.

    # Arguments
    - `linImgs::Vector{LF{Intensity,T,2}}`: A vector of 2D LatticeFields.
    - `idxs::Vector{<:Number}`: Indices corresponding to the LatticeFields in `linImgs`.
    - `roi`: An optional region of interest to consider within each LatticeField.
    - `threshold::Number`: A threshold value for clipping intensities (default: 0.1).

    # Returns
    - A tuple containing the position of the centroid and the orientation angle (theta) of the images.

    The function first calculates the centroid of each LatticeField, either in full or within the specified ROI. 
    It then performs a linear fit to these centroid positions to determine the overall orientation.
    """
function getOrientation(linImgs::Vector{LF{Intensity,T,2}}, idxs::Vector{<:Number}; roi=nothing, threshold::Number=0.1) where {T<:Number}
    if !isnothing(roi)
        linImgs = [i[roi...] for i in linImgs]
    end
    cs = vcat((centroid(i, threshold)' for i in linImgs)...)
    xs, xi, ys, yi = linearFit(idxs, cs[:, 1])..., linearFit(idxs, cs[:, 2])...
    theta = angle(xs + im * ys)
    return [xi, yi], theta
end


end # module ImageProcessing