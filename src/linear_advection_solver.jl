# σ = a dt/dx   =>   dt = dx σ / a
function solve_lin_adv(sys::HyDySys, σ::Float64, a::Float64, t_end::Float64)
  space_order = 2
  t = 0
  N = size(sys.us, 1)
  dt = σ * sys.dx / a
  dx = sys.dx
  println("Started evaluation with following Parameters:")
  println("dt = ", dt)
  println("dx = ", dx)
  println("N = ", N)
  
  # add ghost cells, according to boundary condition
  if sys.bound_cond == :periodic
    ρs = vcat(sys.ρs[end-space_order+1:end], sys.ρs, sys.ρs[1:space_order])
    us = vcat(sys.us[end-space_order+1:end], sys.us, sys.us[1:space_order])
  else
    error("Boundary condition not valid: $(sys.bound_cond)")
    return false
  end

  
  Δρ(ρs, j) = (ρs[j+1] - ρs[j])*(ρs[j] - ρs[j-1]) > 0 ? 
              (ρs[j+1] - ρs[j])*(ρs[j] - ρs[j-1]) / max(ρs[j+1] + ρs[j-1], 0.0001) : 0

  ρ_adv(ρs, us, j, dt_dx) = us[j] > 0 ?
                            ρs[j-1] + 1/2*(2-us[j]*dt_dx) * Δρ(ρs, j-1) :
                            ρs[j]   - 1/2*(2+us[j]*dt_dx) * Δρ(ρs, j)
  
  Fm(ρs, us, j, dt_dx) = ρ_adv(ρs, us, j, dt_dx) * us[j]
  
  while t <= t_end
    println("t = ", t)
    t += dt
    ρs_copy = copy(ρs)
    us_copy = copy(us)

    for j in space_order+1:N+space_order
    #   ρ_adv0 = 0
    #   if us_copy[j] > 0
    #     j_delta = j-1
    #     Δρ = 2 * max(0, (ρs_copy[j_delta+1] - ρs_copy[j_delta]) * (ρs_copy[j_delta] - ρs_copy[j_delta-1]))/(ρs_copy[j_delta+1] - ρs_copy[j_delta-1] + 0.00001)
    #     ρ_adv0 = ρs_copy[j-1] + 0.5+(1-us_copy[j] * dt/dx) * Δρ
    #   else
    #     j_delta = j
    #     Δρ = 2 * max(0, (ρs_copy[j_delta+1] - ρs_copy[j_delta])(ρs_copy[j_delta] - ρs_copy[j_delta-1]))/(ρs_copy[j_delta+1] - ρs_copy[j_delta-1])
    #     ρ_adv0 = ρs_copy[j] - 0.5+(1+us_copy[j] * dt/dx) * Δρ
    #   end
      
    #   ρ_adv1 = 0
    #   if us_copy[j+1] > 0
    #     j_delta = j
    #     Δρ = 2 * max(0, (ρs_copy[j_delta+1] - ρs_copy[j_delta]) * (ρs_copy[j_delta] - ρs_copy[j_delta-1]))/(ρs_copy[j_delta+1] - ρs_copy[j_delta-1] + 0.00001)
    #     ρ_adv0 = ρs_copy[j] + 0.5+(1-us_copy[j+1] * dt/dx) * Δρ
    #   else
    #     j_delta = j+1
    #     Δρ = 2 * max(0, (ρs_copy[j_delta+1] - ρs_copy[j_delta])(ρs_copy[j_delta] - ρs_copy[j_delta-1]))/(ρs_copy[j_delta+1] - ρs_copy[j_delta-1])
    #     ρ_adv0 = ρs_copy[j+1] - 0.5+(1+us_copy[j+1] * dt/dx) * Δρ
    #   end

    #   Fm1 = ρ_adv1 * us_copy[j+1]
    #   Fm0 = ρ_adv0 * us_copy[j]
    #   # ρs[j] = ρs_copy[j] - σ/a * (Fm1 - Fm0)
      ρs[j] = ρs_copy[j] - σ/a * (Fm(ρs_copy, us_copy, j+1, σ/a) - Fm(ρs_copy, us_copy, j, σ/a))
    end
  end

  new_sys = HyDySys(ρs[space_order+1:end-space_order], dx, us[space_order+1:end-space_order], sys.bound_cond)
  return new_sys
end

