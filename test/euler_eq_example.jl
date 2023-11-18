function small_test4()
	#N  = 100
	#dx = 2/(N-1)
	#a  = 1.0
	#σ  = 0.8
	#dt = σ * dx / a


	N = 100
	dx = 1/N
	a = 1.0
	dt = 0.001
	σ = dt * a / dx

	T  = 600.0*dt

	ρs = zeros(N)
	ρs[div(3*N, 8):div(5*N, 8)] = ones(Float64, size(div(3*N, 8):div(5*N, 8), 1))
	us = ones(Float64, N) .* a

	#sys = HyDySys(ρs, dx, us, :periodic, zeros(N), 0.5)
	#new_sys1 = solve_shock_tube(sys, σ, a, T)
	#new_sys2 = solve_lin_adv(sys, σ, a, T)
	#display([sys, new_sys1])
	#display([sys, new_sys2])

	sys_refl = HyDySys(ρs, dx, us, :reflective, zeros(N), 1.4)
	new_sys_refl = solve_euler(sys_refl, σ, a, T)
	display([sys_refl, new_sys_refl])
end


function A2()
	function make_start_sys(x_0, ul, ur)
		N = 100
		dx = 0.01
		γ = 1.4
		dt = 0.001

		σ = 0.99# dt * a / dx
		
		xs = range(0, 1, N)
		ρs = zeros(Float64, N)
		us = zeros(Float64, N)
		ϵs = zeros(Float64, N)
		for i in 1:N
			if xs[i] <= 0.5
				ρs[i], us[i], ϵs[i] = ul
			else
				ρs[i], us[i], ϵs[i] = ur
			end
		end
		return HyDySys(ρs, dx, us, :periodic, ϵs, γ)
	end
	ul = (1.000, 2.5, 0.0)
	ur = (0.125, 2.0, 0.0)
	sys1 = make_start_sys(0.5, ul, ur)
	# println(sys1)
	# solve_lin_adv(sys1, 0.8, 1.0, 0.0000001) |> visualize_system

	frames = [sys1]
	steps = 10
	for i in 1:steps
		append!(frames, [solve_lin_adv(frames[end], 0.999, 1.0, 0.0000001)])
		# append!(frames, [solve_shock_tube(frames[end], 0.999, 1.0, 0.0000001)])
	end
	visualize_system(frames, save_path="media/A2.gif")
	# solve_shock_tube(sys1, 0.8, 0.0, 1.0) |> display
end