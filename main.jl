using Pkg
function installed()
	deps = Pkg.dependencies()
	installs = Dict{String, VersionNumber}()
	for (uuid, dep) in deps
		dep.is_direct_dep || continue
		dep.version === nothing && continue
		installs[dep.name] = dep.version
	end
	return installs
end
# Check if packages are installed, else install them
Packages = ["Plots"]
installed_Packages = keys(installed())
for Package in Packages
	if !(Package in installed_Packages)
		try
			eval(Meta.parse("using $Package"))
		catch
			Pkg.add(Package)
			println("Package $Package was not found. Installation started")
			eval(Meta.parse("using $Package"))
		end
	else
		eval(Meta.parse("using $Package"))
	end
end

include("src/hydrodynamic_system_structs.jl")
include("src/linear_advection_solver.jl")
include("src/shock_tube_solver.jl")
include("src/utilities.jl")
include("test/lin_adv_example.jl")
include("test/euler_eq_example.jl")

# small_test1()
A2()
# Test(Ψ, [-1,1], 1, 4, [40,400],:periodic)