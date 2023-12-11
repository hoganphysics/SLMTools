# This script can be run to add dependencies to the module.
# It should be run from the root directory of the package.
# After doing this, you may have to restart the Julia session to use the new dependencies.
using Pkg
Pkg.activate(".")
# Pkg.add("Images") # add the package
# Pkg.add("Interpolations") # add the package
# Pkg.add("FileIO") # add the package
# Pkg.add("Plots") # add the package
# Pkg.add("FFTW") # add the package
# Pkg.add("Statistics") # add the package
# Pkg.add("FreeTypeAbstraction") # add the package
# Pkg.add("OptimalTransport") # add the package
# Pkg.add("PyCall")
# Pkg.add("ImageMagick")

Pkg.add(["Images", "Interpolations", "FileIO", "Plots", "FFTW", "Statistics",
    "FreeTypeAbstraction", "OptimalTransport", "PyCall", "ImageMagick"])

Pkg.instantiate() # Install the packages
Pkg.precompile() # Precompile the packages
