using Plots

function Base.show(io::IO, ::MIME"text/plain", sys::HyDySys)
	y = sys.ρs
	x = 1:sys.dx:size(y, 1)
	display(plot(x, y))
end
function Base.show(io::IO, ::MIME"text/plain", systems::Vector{HyDySys})
	xs = [[sys.dx .* 1:size(sys.us, 1)] for sys in systems]
	ys = [[sys.ρs] for sys in systems]
	display(xs)
	display(ys)
	display(plot(xs, ys))
end

