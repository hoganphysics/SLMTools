using PyCall


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
function saveAs8BitBMP(image_data::Array{Int64,2}, output_filename::String)
    py_executable_path = PyCall.python
    ENV["PYTHON"] = py_executable_path
    # Read the entire Python script as a string    
    python_script = read("toBbmp.py", String)

    # Execute the Python script
    py"exec($python_script)"

    # Convert Julia array to a format that can be passed to Python
    py_image_data = PyObject(image_data)

    # Call the Python function with the image data and output filename
    py"save_as_bmp($py_image_data, output_filename)"
end
