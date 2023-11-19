using Statistics

function solve_euler(sys::HyDySys, σ::Float64, a::Float64, t_ende::Float64)

	# Identifiziere wichtige Größen
	t       = 0.0
	l       = size(sys.us, 1)
	Δx      = sys.dx
	Δt      = σ * Δx / a
	ρs      = copy(sys.ρs)
	us      = copy(sys.us)
	ϵs      = copy(sys.ϵs)
	γ       = sys.γ
	ordnung = 2 +1

	# Implementiere Formeln zur Lösung nach dem Schema des Opeartor-Splitting


	# *********************************************************
	# ***** 1) Advektionsschritt: A^1 = A^n + Δt*L_1(A^n) *****
	# *********************************************************


	# ******************
	# ***** Dichte *****
	# ******************

	# (2.16)
	function Δρ(ρs, j)
		# Vermeidung doppelter Berechnung
		x = (ρs[j+1]-ρs[j])*(ρs[j]-ρs[j-1])
		divisor = ρs[j+1]-ρs[j-1]
		if x > 0
			if divisor == 0
				@warn "divisor is zero" maxlog=10
				return 2 * x / 0.000001
			elseif divisor != 0
				return 2 * x / (ρs[j+1]-ρs[j-1])
			end
		else
			return 0
		end
	end

	# (2.15)
	function ρ_adv(ρs, us, j, Δt, Δx)
		if us[j] > 0
			return ρs[j-1] + (1/2)*(1-us[j]*Δt/Δx) * Δρ(ρs, j-1)
		else
			return ρs[j]   - (1/2)*(1+us[j]*Δt/Δx) * Δρ(ρs, j)
		end
	end

	# (2.14)
	function Fm(ρs, us, j, Δt, Δx)
		return ρ_adv(ρs, us, j, Δt, Δx) * us[j]
	end

	# (2.13)
	function ρ_new(ρs, us, j, Δt, Δx)
		return ρs[j] - Δt/Δx * (Fm(ρs, us, j+1, Δt, Δx) - Fm(ρs, us, j, Δt, Δx))
	end


	# ******************
	# ***** Impuls *****
	# ******************

	# (2.22)
	function Δu(us, j)
		# Vermeidung doppelter Berechnung
		x = (us[j+1]-us[j])*(us[j]-us[j-1])
		divisor = us[j+1]-us[j-1]
		if x > 0
			if divisor == 0
				@warn "divisor is zero" maxlog=10
				return 2 * x / 0.000001
			elseif divisor != 0
				return 2 * x / divisor
			end
		else
			return 0
		end
	end

	# (2.20)
	function u_adv(us, j, Δt, Δx)
		u_mean = mean([us[j], us[j+1]]) 
		if u_mean > 0 
			return us[j]   + (1/2) * (1 - u_mean * Δt/Δx) * Δu(us, j)
		else
			return us[j+1] - (1/2) * (1 + u_mean * Δt/Δx) * Δu(us, j+1)
		end
	end

	# (2.19)
	function Fl(ρs, us, j, Δt, Δx)
		return mean([Fm(ρs, us, j, Δt, Δx), Fm(ρs, us, j+1, Δt, Δx)]) * u_adv(us, j, Δt, Δx)
	end

	# (2.23)
	function zwischenwert_u(ρs, ρs_new, us, j, Δt, Δx)
		divisor = mean([ρs_new[j-1], ρs_new[j]])
		if divisor == 0
			@warn "divisor is zero" maxlog=10
			return ((us[j]*mean([ρs[j], ρs[j-1]])) - Δt/Δx * (Fl(ρs, us, j, Δt, Δx) - Fl(ρs, us, j-1, Δt, Δx))) / 0.000001
		else
			return ((us[j]*mean([ρs[j], ρs[j-1]])) - Δt/Δx * (Fl(ρs, us, j, Δt, Δx) - Fl(ρs, us, j-1, Δt, Δx))) / divisor
		end
	end


	# *******************
	# ***** Energie *****
	# *******************

	# (2.29)
	function Δϵ(ϵs, j)
		# Vermeidung doppelter Berechnung
		x = (ϵs[j+1]-ϵs[j]) * (ϵs[j]-ϵs[j-1])
		if x > 0
			return 2 * x / (ϵs[j+1]-ϵs[j-1])
		else
			return 0
		end
	end

	# (2.28)
	function ϵ_adv(ϵs, us, j, Δt, Δx) 
		if us[j] > 0
			return ϵs[j-1] + (1/2) * (1 - us[j] * Δt/Δx) * Δϵ(ϵs, j-1)
		else
			return ϵs[j]   - (1/2) * (1 + us[j] * Δt/Δx) * Δϵ(ϵs, j)
		end
	end

	# (2.27)
	function Fe(ρs, us, ϵs, j, Δt, Δx)
		return Fm(ρs, us, j, Δt, Δx) * ϵ_adv(ϵs, us, j, Δt, Δx)
	end

	# (2.30)
	function zwischenwert_ϵ(ρs, ρs_new, us, ϵs, j, Δt, Δx)
		divisor = mean([ρs_new[j-1], ρs_new[j]])
		if divisor == 0
            @warn "divisor is zero" maxlog=10
			return (ϵs[j]*ρs[j] + Δt/Δx * (Fe(ρs, us, ϵs, j+1, Δt, Δx) - Fe(ρs, us, ϵs, j, Δt, Δx))) / 0.000001
		else
			return (ϵs[j]*ρs[j] + Δt/Δx * (Fe(ρs, us, ϵs, j+1, Δt, Δx) - Fe(ρs, us, ϵs, j, Δt, Δx))) / divisor
		end
	end



	# ***************************************************************
	# ***** 2) Kräfte, Druckarbeit: A^{n+1} = A^1 + Δt*L_2(A^1) *****
	# ***************************************************************


	# ******************
	# ***** Impuls *****
	# ******************

	# (2.36)
	function zwischenwert_p(ρs, ϵs, func_γ, j)
		return (func_γ - 1) * ρs[j] * ϵs[j]
	end

	function u_new(ρs_new, us_zwischen, ϵs_zwischen, j ,Δt, Δx, func_γ)
		divisor = mean([ρs_new[j-1], ρs_new[j]])
		if divisor == 0
			@warn "divisor is zero" maxlog=10
			return us_zwischen[j] - Δt * (1/0.00001) * ((func_γ-1)*(ρs_new[j]ϵs_zwischen[j]-ρs_new[j-1]ϵs_zwischen[j-1])/Δx)
		else
			return us_zwischen[j] - Δt * (1/divisor) * ((func_γ-1)*(ρs_new[j]ϵs_zwischen[j]-ρs_new[j-1]ϵs_zwischen[j-1])/Δx)
		end
	end

	# *******************
	# ***** Energie *****
	# *******************   

	function ϵ_new(ρs_new, us_zwischen, ϵs_zwischen, j ,Δt, Δx, func_γ)
		divisor = ρs_new[j+1]
		if divisor == 0
			@warn "divisor is zero" maxlog=10
			return ϵs_zwischen[j] - Δt * (((func_γ-1)*ρs_new[j]*ϵs_zwischen[j])/0.000001) * ((us_zwischen[j+1]-us_zwischen[j])/Δx)
		elseif divisor != 0
			return ϵs_zwischen[j] - Δt * (((func_γ-1)*ρs_new[j]*ϵs_zwischen[j])/ρs_new[j+1]) * ((us_zwischen[j+1]-us_zwischen[j])/Δx)
		end
	end

	function add_ghost_cells!(ρs, us, ϵs, bound_cond, order)
		ρs, us, ϵs = copy.([ρs, us, ϵs])
		if bound_cond == :reflective
			us[1]   = 0.0
			us[end] = 0.0
			# us = vcat(zeros(Float64, order-1), -us[order+2], us, -us[end-1], zeros(Float64, order-1))
			us = vcat(-us[order+1:-1:2], us, us[(end-1):-1:(end-order)])
			ρs = vcat(ρs[order:-1:1], ρs, ρs[(end):-1:(end-order+1)])
			ϵs = vcat(ϵs[order:-1:1], ϵs, ϵs[(end):-1:(end-order+1)])
			return ρs, us, ϵs
		elseif bound_cond == :periodic
			ρs = vcat(ρs[end-order+1:end], ρs, ρs[1:order])
			us = vcat(us[end-order+1:end], us, us[1:order])
			ϵs = vcat(ϵs[end-order+1:end], ϵs, ϵs[1:order])
			return ρs, us, ϵs
		else
			error("Boundary condition not valid: $(bound_cond)")
			return false
		end
	end

	# Beginne nun mit der eigentlichen Berechnung
	while t < t_ende
		# Füge Geister-Zellen basierend auf den Randbedingungen (2.39) ein
		ρs, us, ϵs = add_ghost_cells!(ρs, us, ϵs, sys.bound_cond, ordnung)

		# Fertige Kopien der Systemgrößen an, damit bei der Berechnung nichts durcheinander kommt
		ρs_copy, us_copy, ϵs_copy = copy(ρs), copy(us), copy(ϵs)


		# Berechnung des Advektionsschritts, die updates erfolgen nach (2.31):
		# ρ → ρ^{n+1}
		# u → u^1
		# ϵ → ϵ^1
		# Berechne zuerst die neue Dichte des Systems:
		for j in (ordnung+1):(l + ordnung - 1)
			ρs[j] = ρ_new(ρs_copy, us_copy, j, Δt, Δx)
		end

		# Berechne nun die Zwischenwerte der Geschwindigkeit und Energie
		for j in (ordnung+1):(l + ordnung - 1)
			# If-Schleife damit [j=3...N+1] eingehalten wird
			if (j >= (ordnung+2))
				us[j] = zwischenwert_u(ρs_copy, ρs, us_copy, j, Δt, Δx)
			end
			ϵs[j] = zwischenwert_ϵ(ρs_copy, ρs, us_copy, ϵs_copy, j, Δt, Δx) 
		end

		#println("Dichte zu t=", t, ": " , ρs)
		#println("Zwischen-Geschwindigkeit:", us)
		#println("Zwischen-Energie:", ϵs)

		# Füge Geister-Zellen basierend auf den Randbedingungen (2.39) ein
		ρs, us, ϵs = ρs[ordnung+1:end-ordnung], us[ordnung+1:end-ordnung], ϵs[ordnung+1:end-ordnung]
		ρs, us, ϵs = add_ghost_cells!(ρs, us, ϵs, sys.bound_cond, ordnung)


		# Berechnung der Kräfte und Druckarbeit, hier werden u und ϵ geupdated
		# u → u^1 → u^{n+1}
		# ϵ → ϵ^1 → ϵ^{n+1}
		for j in (ordnung + 1):(l + ordnung - 1)
			# If-Schleife damit [j=3...N+1] eingehalten wird
			if j >= (ordnung+2)
				us[j] = u_new(ρs, us, ϵs, j, Δt, Δx, γ)
			end
			ϵs[j] = ϵ_new(ρs, us, ϵs, j, Δt, Δx, γ)
		end

		#println("Update-Geschwindigkeit:", us)
		#println("Update-Energie:", ϵs)

		# Entferne Giisterzellen
		ρs, us, ϵs = ρs[ordnung+1:end-ordnung], us[ordnung+1:end-ordnung], ϵs[ordnung+1:end-ordnung]

		# Update Zeitschritt
		t += Δt
	end
	return HyDySys(ρs, Δx, us, sys.bound_cond, ϵs, γ)
end