# This script can be run to add dependencies to the module.
# It should be run from the root directory of the package. although currently it does that by setting the base_dir variable.
# After doing this, you may have to restart the Julia session to use the new dependencies.
using Pkg
base_dir = "Z:\\home\\Vandy\\code\\julia\\SLMTools"
cd(base_dir)
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
# "StochasticOptimalTransport"
Pkg.add(["Images", "Interpolations", "FileIO", "Plots", "FFTW", "Statistics",
    "FreeTypeAbstraction", "OptimalTransport", "ImageMagick","Revise"])

Pkg.instantiate() # Install the packages
Pkg.precompile() # Precompile the packages

# base_dir = "Z:\\home\\Vandy\\code\\julia\\SLMTools"
# cd(base_dir)
# Pkg.activate(".")
# Pkg.add("SLMTools") # add the package