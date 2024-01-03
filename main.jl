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
Packages = ["Plots" "LaTeXStrings" "Statistics"]
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
include("src/utilities.jl")
include("src/lin_adv_solver.jl")
include("test/lin_adv_example.jl")
include("src/euler_solver.jl")



# Löse A1
solve_adv()

# Löse A2
solve_euler(100, 1.4)