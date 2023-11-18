function small_test4()
	N  = 100
	dx = 2/(N-1)
	a  = 1.0
	σ  = 0.8
	dt = σ * dx / a

	#N = 100
	#dx = 1/N
	#a = 1.0
	#dt = 0.001
	#σ = dt * a / dx

	T  = 70.0*dt

	ρs = zeros(N)
	ρs[div(3*N, 8):div(5*N, 8)] = ones(Float64, size(div(3*N, 8):div(5*N, 8), 1))
	us = ones(Float64, N) .* a

	sys = HyDySys(ρs, dx, us, :periodic, zeros(N), 0.5)
	new_sys1 = solve_shock_tube(sys, σ, a, T)
	new_sys2 = solve_lin_adv(sys, σ, a, T)
	display([sys, new_sys1])
	display([sys, new_sys2])

	#sys_refl = HyDySys(ρs, dx, us, :reflective, zeros(N), 1.4)
	#new_sys_refl = solve_shock_tube(sys_refl, σ, a, T)
	#display([sys_refl, new_sys_refl])
end