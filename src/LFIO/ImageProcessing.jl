#=
Functions for working with real images in the context of beam estimation.  There are many conventions we have adopted in the interest of 
    developing a convenient workflow for this generally challenging and annoying problem. 

Requires:
    using FileIO: load, save
    using Images: Gray, RGB, Colorant
=#

export getImagesAndFilenames, imageToFloatArray, itfa, castImage, loadDir, parseFileName, parseStringToNum, getOrientation, dualate, linearFit, savePhase, saveBeam, saveAs8BitBMP

"""
    castImage(T::DataType, img::Matrix{<:Colorant}, L::Lattice{2}, flambda::Real)

    Cast an image of colorants into a specified data type and create a LatticeField representation.

    This function takes an image matrix `img`, which consists of colorants (such as RGB values), and casts it 
    into a specific data type `T`. It then creates and returns a `LatticeField` object of type `T` based on the 
    transformed image, the specified lattice `L`, and a wavelength parameter `flambda`.

    # Arguments
    - `T::DataType`: The FieldVal type to which the image should be cast.
    - `img::Matrix{<:Colorant}`: A matrix representing the image, where each element is a colorant (e.g., an RGB value).
    - `L::Lattice{2}`: A 2-dimensional lattice that defines the spatial properties of the resulting LatticeField.
    - `flambda::Real`: A real number representing the wavelength parameter for the LatticeField.

    # Returns
    - `LF{T}`: A LatticeField of type `T`, representing the transformed image.

    # Example
    ```julia
    julia> img = load("path/to/image.png") # Assume this returns a 100x100 Matrix{RGB}
    julia> L = Lattice((100, 100)) # Define a 2D lattice
    julia> flambda = 1.064e5 # Scale factor for a 10 cm lens and 1064 nm wavelength. 
    julia> transformed_img = castImage(Intensity, img, L, flambda)
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
    - `T::DataType`: The FieldVal type to assign the imported images (default: `Intensity`).
    - `outType::DataType`: The numerical data type of the output LatticeField arrays (default: `Float64`).
    - `L::Union{Lattice{2},Nothing,Real,Tuple{Real,Real}}`: Lattice defining spatial properties. If `Real` or `Tuple`, it's interpreted as pixel spacing.
    - `flambda::Real`: Wavelength parameter for the LatticeFields (default: 1.0).
    - `cue::Union{String,Char,Nothing}`: Optional string or character used as a cue to parse filenames (default: `nothing`).  If `nothing`, then the entire filename (except the extension) is parsed. 
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
    # Converts a numeric string to a number.  Numeric strings can only contain numerals, commas, and minus signs. 
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
    # The numeric string can only contain numerals, a minus sign, and a comma, which stands in for a decimal point.
    # By default, if there's a comma the numeric string is converted to a float, and otherwise it is converted
    # to an integer.  
    look in [:before, :b, "before", "b", 'b', :after, :a, "after", "a", 'a'] || throw(ValueError("Unrecognized look value."))
    look in [:before, :b, "before", "b", 'b'] ? (look = :b) : (look = :a)
    ext = findlast('.', name)
    name = '.' * name[1:ext-1] * '.'    # The periods are to simplify finding the beginning/end of the numeric string. 
    idx = findlast(cue, name)
    if look == :a
        j = findfirst((x -> !(x in [" 0123456789,-"...])), name[idx[end]+1:end]) + idx[end]
        s = name[idx[end]+1:j-1]
        return parseStringToNum(s, outType=outType)
    else
        j = findlast((x -> !(x in [" 0123456789,-"...])), name[1:idx[1]-1])
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
	linearFit(xs::Vector{<:Number},ys::Vector{<:Number})
	
	Returns (slope,intercept) of the least squares regression linear fit to the points with x-coordinates `xs` and
	y coordinates `ys`. 
	"""
function linearFit(xs::Vector{<:Number},ys::Vector{<:Number})
    # Fits a line to the points (x,y) with x in xs, y in ys, via linear regression.  Returns (slope, intercept).
    ((hcat(xs, ones(length(xs))) \ ys)...,)
end

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

"""
    dualate(f::LF{S,T,2}, L::Lattice{2}, center::Vector{<:Real}, θ::Real, flambda::Real=1.0; 
            roi=nothing, interpolation=cubic_spline_interpolation, naturalize=false, bc=zero(T)) 
            where {S<:FieldVal, T<:Number}

    Corrects for offset and tilt in a LatticeField and interpolates it onto a lattice dual to `L`.

    This function adjusts the LatticeField `f` for any offset and tilt based on `center` and `θ`, and then performs interpolation 
    onto a new lattice that is dual to `L`. The dual lattice is calculated considering the scaling factor `flambda`.

    # Arguments
    - `f::LF{S,T,2}`: A two-dimensional LatticeField to be dualated.
    - `L::Lattice{2}`: The target lattice for dualation.
    - `center::Vector{<:Real}`: The center around which the dualation is performed.
    - `θ::Real`: The angle of rotation for tilt correction.
    - `flambda::Real`: Scaling factor (default: 1.0).
    - `roi`: Optional region of interest to be considered within `f`.
    - `interpolation`: Interpolation method to be used (default: cubic spline interpolation).
    - `naturalize`: Boolean flag to naturalize the lattice (default: false).
    - `bc`: Boundary condition for extrapolation (default: `zero(T)`).

    # Returns
    - A new `LF{S}` representing the dualated LatticeField.
    """
function dualate(f::LF{S,T,2},L::Lattice{2},center::Vector{<:Real},θ::Real,flambda::Real=1.0; 
        roi=nothing, interpolation=cubic_spline_interpolation, naturalize=false,bc=zero(T)) where {S<:FieldVal,T<:Number}
    # Corrects for offset and tilt of f.L and then interpolates f onto a lattice dual to L. 
    length(center) == length(f.L) || error("Incompatible lengths for center and f.L.")
    if !isnothing(roi)
        f = f[roi...]
    end
    L0 = Tuple(f.L[i] .- center[i] for i=1:2)
    interp = interpolation(L0, f.data, extrapolation_bc=bc)
    dL = dualShiftLattice(L,flambda)
    Lx0, Lxs, Ly0, Lys = dL[1][1], step(dL[1]), dL[2][1], step(dL[2])
    R = [cos(θ) -sin(θ) ; sin(θ) cos(θ)]
    r0, dx, dy = (R * [Lx0, Ly0]), Lxs * R[:,1], Lys * R[:,2]
    if naturalize
        return LF{S}([interp((r0 + dx*I[1] + dy*I[2])...) for I in CartesianIndices(length.(dL))],natlat(length.(dL)),1.0);
    else
        return LF{S}([interp((r0 + dx*I[1] + dy*I[2])...) for I in CartesianIndices(length.(dL))],dL,flambda);
    end
end
dualate(fs::Vector{<:LF{S}},L::Lattice{2},center::Vector{<:Real},theta::Real,flambda::Real=1.0;kwargs...) where {S<:FieldVal} = 
    [dualate(f,L,center,theta,flambda;kwargs...) for f in fs]


"""
    savePhase(A::AbstractArray{ComplexF64,2}, name::String)

    Save the phase of a complex array as a grayscale image.

    # Arguments
    - `A::AbstractArray{ComplexF64,2}`: A 2D array of complex numbers whose phase is to be saved.
    - `name::String`: The filename for the saved image.

    This function normalizes the phase of the complex array `A` to the range [0, 1] by adding π, dividing by 2π, and then converting it to grayscale before saving the image with the given `name`.
"""
function savePhase(A::AbstractArray{ComplexF64,2}, name::String)
    save(name, Gray.((angle.(A) .+ pi) ./ (2pi)))
end

"""
    savePhase(x::LF{ComplexPhase}, name::String)

    Save the phase information of a LatticeField with complex phase values to a file.

    This function extracts the phase data from a `LF{ComplexPhase}` object and saves it using the specified file name.

    # Arguments
    - `x::LF{ComplexPhase}`: A LatticeField object containing complex phase values.
    - `name::String`: The name of the file to save the phase data.

    The function delegates the saving process to another `savePhase` method tailored for handling the raw data of `x`.
    """
function savePhase(x::LF{ComplexPhase}, name::String)
    savePhase(x.data, name)
end

"""
    savePhase(x::LF{RealPhase}, name::String)

    Save the phase information of a LatticeField with real phase values as an image file.

    This function takes a `LF{RealPhase}` object, normalizes the phase data to the range [0,1), and saves it as an image 
    using the specified file name. The normalization is performed by taking the modulus of the phase data with 1.

    # Arguments
    - `x::LF{RealPhase}`: A LatticeField object containing real phase values.
    - `name::String`: The name of the image file to save the normalized phase data.

    The saved image visualizes the phase information by mapping phase values to grayscale intensities.
    """
function savePhase(x::LF{RealPhase}, name::String)
    save(name, Gray.(mod.(x.data, 1)))
end

"""
    saveBeam(beam::Matrix{ComplexF64}, name::String, data=[:beamCsv, :angleCsv, :anglePng, :negativeAnglePng]; dir=pwd())

    Save various representations of a complex beam array to files.

    # Arguments
    - `beam::Matrix{ComplexF64}`: A matrix representing a complex beam.
    - `name::String`: The base filename for the saved files.
    - `data`: An array of symbols indicating which types of data to save. Options are `:beamCsv`, `:angleCsv`, `:anglePng`, and `:negativeAnglePng`.
    - `dir`: The directory in which to save the files, defaults to the present working directory.

    # Details
    This function can save the complex beam and its angle as CSV files and PNG images. It appends the current date to the filenames and allows for saving the negative of the phase as well.
    """
function saveBeam(beam::Matrix{ComplexF64}, name::String, data=[:beamCsv, :angleCsv, :anglePng, :negativeAnglePng]; dir=pwd())
    if :beamCsv in data
        writedlm(dir * "\\" * name * "-Beam-" * string(now())[1:10] * ".csv", beam, ',')
    end
    if :angleCsv in data
        writedlm(dir * "\\" * name * "-Angle-" * string(now())[1:10] * ".csv", angle.(beam), ',')
    end
    if :anglePng in data
        savePhase(beam, dir * "\\" * name * "-Angle-" * string(now())[1:10] * ".png")
    end
    if :negativeAnglePng in data
        savePhase(beam, dir * "\\" * name * "-NegativeAngle-" * string(now())[1:10] * ".png")
    end
end


"""
    saveAs8BitBMP(image_data::Array{Int64,2}, output_filename::String)

    Execute a Python script and save an image from Julia array data.

    # Arguments
    - `image_data::Array{Int64,2}`: A 2D array of integers representing the image data.
    - `output_filename::String`: The filename for the saved image.

    This function executes python code which saves an image from the given Julia array data. The Python script is read from the file `script` and the image data is passed to the Python function `save_as_bmp` along with the output filename.
    """
# function saveAs8BitBMP(image_data::Array{Int64,2}, output_filename::String)
#     py_executable_path = PyCall.python
#     ENV["PYTHON"] = py_executable_path
#     # Read the entire Python script as a string    
#     python_script = read("to8Bbmp.py", String)

#     # Execute the Python script
#     py"exec($python_script)"

#     # Convert Julia array to a format that can be passed to Python
#     py_image_data = PyObject(image_data)

#     # Call the Python function with the image data and output filename
#     py"save_as_bmp($py_image_data, $output_filename)"
# end
