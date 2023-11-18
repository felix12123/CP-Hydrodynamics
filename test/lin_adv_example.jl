using Plots

function small_test1()
	N  = 40
	dx = 2/(N-1)
	a  = 1.0
	σ  = 0.8

	dt = σ * dx / a
	println("dt = ", dt)
	T  = 2.0

	ρs = zeros(Float64, N)
	ρs[div(3*N, 8):div(5*N, 8)] = ones(Float64, size(div(3*N, 8):div(5*N, 8), 1))
	us = ones(Float64, N) .* a
	
	sys = HyDySys(ρs, dx, us, :periodic)
	new_sys1 = solve_lin_adv(sys, σ, a, T)
	
	steps = T/dt
	frames = [sys]
	for i in 1:steps
		append!(frames, [solve_lin_adv(frames[end], σ, a, dt)])
	end

	visualize_system(frames, disp=false, save_path="media/A1_1.gif")
end

function small_test2()
	N  = 2<<8
	dx = 2/(N-1)
	a  = 1.0
	σ  = 0.8

	dt = σ * dx / a
	T  = 2.0  #dx*N/a / 500

	gauss(x, μ, σ) = 1/(sqrt(2pi) * σ) * exp(-1/2*(x-μ)^2/σ^2)

	xs = (1:N) .* dx
	ρs = gauss.(xs, N/2*dx, N/25*dx)
	us = ones(Float64, N) .* a

	sys = HyDySys(ρs, dx, us, :periodic)
	new_sys = solve_lin_adv(sys, σ, a, T)
	new_sys1 = solve_shock_tube(sys, σ, a, T)

	# display([sys, new_sys])
	display([sys, new_sys])
	display([sys, new_sys1])
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
		dx     = (interval[2]-interval[1])/(gs-1)
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

  
		#plot(Gitter, Gitter_analytisch, title="Vorhersage für Gridsize ", linewidth=3, gs, label="", dp=300, color= :black)
		#display(plot!(Gitter, Gitter_HyDySys, label="", dp=300, color=:red))
	end
end