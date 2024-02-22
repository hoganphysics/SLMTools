using Pkg
# what is this for? 
base_dir ="Z:\\home\\Vandy\\code\\julia"
cd(base_dir)
# Package name
mypkgname = "SLMTools"
full_path = joinpath(base_dir, mypkgname)

cd(base_dir)
Pkg.activate()
Pkg.develop(path=full_path)