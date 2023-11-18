using Plots

function Base.show(io::IO, ::MIME"text/plain", sys::HyDySys)
	y = sys.ρs
	x = (0:size(y, 1)-1) .* sys.dx |> collect
	display(plot(x, y))
end
function Base.show(io::IO, ::MIME"text/plain", systems::Vector{HyDySys})
	xs = [[sys.dx .* (0:size(sys.us, 1)-1)  |> collect] for sys in systems]
	ys = [[sys.ρs] for sys in systems]
	display(plot(xs, ys, alpha=0.75))
end


function visualize_system(sys::HyDySys; variable=:ρ, title="", disp=true, save_path="")
	xs = eachindex(sys.ρs) .* sys.dx
	if variable == :ρ
		ys = sys.ρs
	elseif variable == :ϵ
		ys = sys.ϵs
	elseif variable == :u
		ys = sys.us
	else
		error("The selected quantity $variable is not valid.")
	end

	plt = plot(xs, ys, title=title, xlabel="position", ylabel=string(variable), label="")
	if save_path != ""
		if ispath(save_path)
			savefig(plt, save_path)
		else
			@warn "path $save_path is not a valid path. picture not saved."
		end
	end

	if disp
		display(plt)
	else
		return plt
	end
end

function visualize_system(sys::Vector{HyDySys}; variable=:ρ, title="", disp=true, save_path="", duration=5)
	pics = Vector{Plots.Plot}(undef, size(sys, 1))
	anim = @animate for i in eachindex(pics)
		pics[i] = visualize_system(sys[i], variable=variable, title=title, disp=false, save_path="")
	end

	if save_path != ""
		gif1 = gif(anim, save_path, fps=size(sys, 1)/5)
		if disp
			display(gif1)
		end
	end
			
	if disp
		display(anim)
	end

end

