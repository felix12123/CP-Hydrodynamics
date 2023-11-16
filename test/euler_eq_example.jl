function small_test4()
	N  = 100
	dx = 2/(N-1)
	a  = 1.0
	σ  = 0.8

	dt = σ * dx / a
	T  = 5.0*dt

	ρs = zeros(Float64, N)
	ρs[div(3*N, 8):div(5*N, 8)] = ones(Float64, size(div(3*N, 8):div(5*N, 8), 1))
	us = ones(Float64, N) .* a

	sys = HyDySys(ρs, dx, us, :periodic, zeros(Float64, N), 0.1)
	new_sys1 = solve_shock_tube(sys, σ, a, T)

	display([sys, new_sys1])
end
