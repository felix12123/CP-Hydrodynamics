function small_test1()
	N  = 40
	dx = 2/N
	a  = 1.0
	σ  = 0.00001

	dt = σ * dx / a
	T  = 4.0

	pack_size = div(N, 4)
	# ρs = vcat(zeros(Float64, floor(Int, (N-pack_size)/2)), ones(Float64, pack_size), zeros(Float64, ceil(Int, (N-pack_size)/2)))
	ρs = zeros(Float64, N)
	ρs[div(3*N, 8):div(5*N, 8)] = ones(Float64, size(div(3*N, 8):div(5*N, 8), 1))
	us = ones(Float64, N) .* a * 2

	sys = HyDySys(ρs, dx, us, :periodic)
	new_sys = solve_lin_adv(sys, σ, a, T)

	display([sys, new_sys])
	# display(new_sys)
end

function small_test2()
	N  = 2<<8
	dx = 100/N
	a  = 1.0
	σ  = 0.008

	dt = σ * dx / a
	T  = dt * 1000 #dx*N/a / 500

	gauss(x, μ, σ) = 1/(sqrt(2pi) * σ) * exp(-1/2*(x-μ)^2/σ^2)

	xs = (1:N) .* dx
	ρs = gauss.(xs, 50, 5)
	us = ones(Float64, N) .* a

	sys = HyDySys(ρs, dx, us, :periodic)
	new_sys = solve_lin_adv(sys, σ, a, T)

	# display([sys, new_sys])
	display([sys, new_sys])
end

# Testing:

# Implementation der Funktoin Ψ(x,t)
function Ψ(x, a, t, interval)
	if abs(mod(x - a * t, (interval[2]-interval[1])/2)) <= 1/3
		return 1
	elseif 1/3 < abs(mod(x - a * t, (interval[2]-interval[1])/2)) <= 1
		return 0
	end
end
  
# Check für das σ = 0.8, Ψ(x, t=4) und 40 Gitterpunkte im Intervall [-1,1]
function Test(Psi_func, interval, a, t, gridsize, bound_cond)
  
	# Führe Test für alle angegebenen Gittergrößen durch
	for gs in gridsize
  
		# Generiere Gitter
		dx     = (interval[2]-interval[1])/gs
		Gitter = interval[1] : dx : interval[2] |> collect
		Gitter_analytisch = copy(Gitter)
		
		# Generiere für den Test das Array der Geschwindigkeiten
		l = length(Gitter)
		us = fill(1.0, l)
		ρs = Vector{Float64}(undef, l)
  
		# Berechne analyt. Lösung und Startarray von ρs
		for i in eachindex(Gitter)
		  # Berechne analyt. Lösung für t=4
		  Gitter_analytisch[i] = Psi_func(Gitter[i], a, t, interval)
		  # Berechne ρs zu t=0
		  ρs[i]                = float(Psi_func(Gitter[i], a, 0, interval))
		end

		#display(HyDySys(ρs, dx, us, bound_cond))
		Gitter_HyDySys = solve_lin_adv(HyDySys(ρs, dx, us, bound_cond), 0.8, 1.0, 4.0)
		display(Gitter_HyDySys)

		if gs == 40
			println(Gitter_HyDySys.ρs)
		end
  
  
		#plot(Gitter, Gitter_analytisch, title="Vorhersage für Gridsize ", linewidth=3, gs, label="", dp=300, color= :black)
		#display(plot!(Gitter, Gitter_HyDySys, label="", dp=300, color=:red))
	end
end