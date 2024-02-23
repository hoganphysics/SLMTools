using Pkg
# this should be run after the addDependencies.jl script
base_dir ="Z:\\home\\Vandy\\code\\julia"
cd(base_dir)
# Package name
mypkgname = "SLMTools"
full_path = joinpath(base_dir, mypkgname)

Pkg.activate()
Pkg.develop(path=full_path)

# Template(interactive=true)("PackageA")