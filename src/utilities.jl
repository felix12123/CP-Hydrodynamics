using Plots

function Base.show(io::IO, ::MIME"text/plain", sys::HyDySys)
	y = sys.ρs
	x = (0:size(y, 1)-1) .* sys.dx |> collect
	display(plot(x, y))
end
function Base.show(io::IO, ::MIME"text/plain", systems::Vector{HyDySys})
	xs = [[sys.dx .* (0:size(sys.us, 1)-1)  |> collect] for sys in systems]
	ys = [[sys.ρs] for sys in systems]
	display(plot(xs, ys))
end

