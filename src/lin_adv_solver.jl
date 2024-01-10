function solve_advection(sys::HyDySys, σ::Float64, a::Float64, t_end::Float64)

	order::Int	 = 2
	t::Float64   = 0.0
	N::Int       = size(sys.ρs, 1)
	Δx::Float64  = sys.dx
	Δt::Float64  = (σ * Δx) / a
	
	ρs::Vector{Float64} = copy(sys.ρs)
	us::Vector{Float64} = copy(sys.us)

	function Δρ(ρs, j)																				# (2.16)
		numer = (ρs[j+1]-ρs[j])*(ρs[j]-ρs[j-1])
		denom = ρs[j+1]-ρs[j-1]
		if numer > 0
			if denom == 0
				return 2 * numer / 0.01
			else
				return 2 * numer / denom
			end
		else
			return 0
		end
	end

	function ρ_adv(ρs, us, j, Δt, Δx)																# (2.15)
		if us[j] > 0
			return ρs[j-1] + (1/2)*(1-us[j]*Δt/Δx) * Δρ(ρs, j-1)
		else
			return ρs[j]   - (1/2)*(1+us[j]*Δt/Δx) * Δρ(ρs, j)
		end
	end

	function Fm(ρs, us, j, Δt, Δx)																	# (2.14)
		return ρ_adv(ρs, us, j, Δt, Δx) * us[j]
	end

	function ρ_new(ρs, us, j, Δt, Δx)																# (2.13)
		return ρs[j] - Δt/Δx * (Fm(ρs, us, j+1, Δt, Δx) - Fm(ρs, us, j, Δt, Δx))
	end

	function add_ghost_cells!(ρs, us, bound_cond, order)
		ρs, us = deepcopy.([ρs, us])
		if bound_cond == :reflective
			us[1]   = 0.0
			us[end] = 0.0
			# us = vcat(zeros(Float64, order-1), -us[order+2], us, -us[end-1], zeros(Float64, order-1))
			us = vcat(-us[order+1:-1:2], us, -us[(end-1):-1:(end-order)])
			ρs = vcat(ρs[order:-1:1], ρs, ρs[(end):-1:(end-order+1)])
			return ρs, us
		elseif bound_cond == :periodic
			ρs = vcat(ρs[end-order+1:end], ρs, ρs[1:order])
			us = vcat(us[end-order+1:end], us, us[1:order])
			return ρs, us
		else
			error("Boundary condition not valid: $(bound_cond)")
			return false
		end
	end

	while t < t_end
		ρs, us = add_ghost_cells!(ρs, us, sys.bound_cond, order)

		ρs_copy = copy(ρs)
		us_copy = copy(us)

		for j in (order+1):(N+order)
			# ρs[j] = ρs_copy[j] - dt/dx * (Fm(ρs_copy, us_copy, j+1, dt, dx) - Fm(ρs_copy, us_copy, j, dt, dx))
			ρs[j] = ρ_new(ρs_copy, us_copy, j, Δt, Δx)
		end
		t += Δt

		ρs = ρs[order+1:end-order]
		us = us[order+1:end-order]
	end
	new_sys = HyDySys(ρs, Δx, us, sys.bound_cond)
	return new_sys
end