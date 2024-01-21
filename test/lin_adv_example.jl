using BenchmarkTools

function solve_adv()

	function make_start_sys(xs, N, Δx, γ)															# Erstelle Startsystem
	
		ρs = zeros(Float64, N)																		# Initialisiere Systemgrößen
		ϵs = ones(Float64, N)
		ps = ones(Float64, N)
		us = ones(Float64, N+1) .* -1

		for j in 1:N																				# Startkonfiguraiton für ρs
			if abs(xs[j]) <= 1/3
				ρs[j] = 1.0
			else
				ρs[j] = 0.0
			end
		end
		
		return HyDySys(ρs, Δx, us, :periodic, ϵs, 1.4, ps)											# Gebe System aus
	end

	N40		 = 40																					# Parameter aus der Aufgabenstellung ↓
	N400	 = 400
	σ  		 = 0.8
	a  		 = 1.0
	Δx40	 = 2/(N40-1)
	Δt40	 = (σ*Δx40)/a
	Δx400    = 2/(N400-1)
	Δt400    = (σ*Δx400)/a
	t_start  = 0.0
	t_end4   = 4.0
	t_end400 = 40.0
	xs40     = collect(range(-1, 1, N40))
	xs400    = collect(range(-1, 1, N400))
	
	start_sys40 = make_start_sys(xs40, N40, Δx40, 1.4)												# Erste Startkonfiguration
	start_sys400 = make_start_sys(xs400, N400, Δx400, 1.4)

	ρs_ana40 = copy(start_sys40.ρs)																	# Analytische Lsg. für ρs nach t_end Zeit (Nach C. Schäfer in beiden Fällen gleich)
	for j in 1:N40
		if abs(xs40[j] - 4 * Δx40) <= 1/3
			ρs_ana40[j] = 1.0
		else
			ρs_ana40[j] = 0.0
		end
	end

	ρs_ana400 = copy(start_sys400.ρs)																# Analytische Lsg. für ρs nach t_end Zeit (Nach C. Schäfer in beiden Fällen gleich)
	for j in 1:N400
		if abs(xs400[j] - 4 * Δx40) <= 1/3
			ρs_ana400[j] = 1.0
		else
			ρs_ana400[j] = 0.0
		end
	end

	# N=40, t=0
	evolved_sys_N40_t0    = solve_advection(start_sys40, σ, a, t_start)
	# N=40, t=4
	@time evolved_sys_N40_t4    = solve_advection(start_sys40, σ, a, t_end4)
	# N=400, t=0
	evolved_sys_N400_t0   = solve_advection(start_sys400, σ, a, t_start)
	# N=400, t=10
	@time evolved_sys_N400_t400 = solve_advection(start_sys400, σ, a, t_end400)

	# Analytischer vs. numerischer Plot zur Zeit N = 40 & t = 0
	# plot(xs40, start_sys40.ρs, label="Analytische Lsg.", linewidth=2, linealpha=0.4, linecolor = :darkblue, dpi= 300, title="", xlabel="x", ylabel=L"\rho", background_color_legend = nothing, fg_legend = :transparent, grid=false)
	plot(xs40, start_sys40.ρs, label="Analytische Lsg.", linewidth=2, linealpha=0.4, linecolor = :darkblue, dpi= 300, title="", xlabel="x", ylabel=L"\rho")
	scatter!(xs40, evolved_sys_N40_t0.ρs, label="Numerische Lsg.", marker=:xcross, markersize=2, markerstrokewidth=2, markercoloer= :orange)
	plot!(xs40, evolved_sys_N40_t0.ρs, label="", linewidth=1, linealpha=0.2, linecolor = :black)
	savefig("media/A1_N40_t0.png")

	# Analytischer vs. numerischer Plot zu N = 40 & t = 4
	# plot(xs40, ρs_ana40, label="Analytische Lsg.", linewidth=2, linealpha=0.4, linecolor = :darkblue, dpi= 300, title="", xlabel="x", ylabel=L"\rho", background_color_legend = nothing, fg_legend = :transparent, grid=false)
	plot(xs40, ρs_ana40, label="Analytische Lsg.", linewidth=2, linealpha=0.4, linecolor = :darkblue, dpi= 300, title="", xlabel="x", ylabel=L"\rho")
	scatter!(xs40, evolved_sys_N40_t4.ρs, label="Numerische Lsg.", marker=:xcross, markersize=2, markerstrokewidth=2, markercoloer= :orange)
	plot!(xs40, evolved_sys_N40_t4.ρs, label="", linewidth=1, linealpha=0.2, linecolor = :black)
	savefig("media/A1_N40_t4.png")

	# Analytischer vs. numerischer Plot zu N = 400 & t = 0
	# plot(xs400, start_sys400.ρs, label="Analytische Lsg.", linewidth=2, linealpha=0.4, linecolor = :darkblue, dpi= 300, title="", xlabel="x", ylabel=L"\rho", background_color_legend = nothing, fg_legend = :transparent, grid=false)
	plot(xs400, start_sys400.ρs, label="Analytische Lsg.", linewidth=2, linealpha=0.4, linecolor = :darkblue, dpi= 300, title="", xlabel="x", ylabel=L"\rho")
	scatter!(xs400, evolved_sys_N400_t0.ρs, label="Numerische Lsg.", marker=:xcross, markersize=2, markerstrokewidth=2, markercoloer= :orange)
	plot!(xs400, evolved_sys_N400_t0.ρs, label="", linewidth=1, linealpha=0.2, linecolor = :black)
	savefig("media/A1_N400_t0.png")

	# Analytischer vs. numerischer Plot zu N = 40 & t = 4
	# plot(xs400, ρs_ana400, label="Analytische Lsg.", linewidth=2, linealpha=0.4, linecolor = :darkblue, dpi= 300, title="", xlabel="x", ylabel=L"\rho", background_color_legend = nothing, fg_legend = :transparent, grid=false)
	plot(xs400, ρs_ana400, label="Analytische Lsg.", linewidth=2, linealpha=0.4, linecolor = :darkblue, dpi= 300, title="", xlabel="x", ylabel=L"\rho")
	scatter!(xs400, evolved_sys_N400_t400.ρs, label="Numerische Lsg.", marker=:xcross, markersize=2, markerstrokewidth=2, markercoloer= :orange)
	plot!(xs400, evolved_sys_N400_t400.ρs, label="", linewidth=1, linealpha=0.2, linecolor = :black)
	savefig("media/A1_N400_t400.png")
end