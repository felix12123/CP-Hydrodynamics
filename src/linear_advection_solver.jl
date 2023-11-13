# σ = a dt/dx   =>   dt = dx σ / a
function solve_lin_adv(sys::HyDySys, σ::Float64, a::Float64, t_end::Float64)
  space_order = 2
  t = 0
  N = size(sys.us, 1)
  dt = σ * sys.dx / a
  dx = sys.dx

  
  # add ghost cells, according to boundary condition
  if sys.bound_cond == :periodic
    ρs = vcat(sys.ρs[end-space_order+1:end], sys.ρs, sys.ρs[1:space_order])
    us = vcat(sys.us[end-space_order+1:end], sys.us, sys.us[1:space_order])
  else
    error("Boundary condition not valid: $(sys.bound_cond)")
    return false
  end

  
  # Define Functions that are needed for better readability of the time step update
  function Δρ(ρs, j)
    if (ρs[j+1] - ρs[j])*(ρs[j] - ρs[j-1]) > 0
      ( ρs[j+1] - ρs[j])*(ρs[j] - ρs[j-1]) / max(ρs[j+1] - ρs[j-1], 0.0001)
    else
      0
    end
  end

  function ρ_adv(ρs, us, j, dt_dx)
    if us[j] > 0
      return ρs[j-1] + 1/2*(1-us[j]*dt_dx) * Δρ(ρs, j-1)
    else
      return ρs[j]   - 1/2*(1+us[j]*dt_dx) * Δρ(ρs, j)
    end
  end
  
  Fm(ρs, us, j, dt_dx) = ρ_adv(ρs, us, j, dt_dx) * us[j]
  
  # Start the actual simulation
  while t <= t_end
    ρs_copy = copy(ρs)
    us_copy = copy(us)
    for j in space_order+1:N+space_order
      ρs[j] = ρs_copy[j] - dt/dx * (Fm(ρs_copy, us_copy, j+1, dt/dx) - Fm(ρs_copy, us_copy, j, dt/dx))
    end
    t += dt
  end

  new_sys = HyDySys(ρs[space_order+1:end-space_order], dx, us[space_order+1:end-space_order], sys.bound_cond)
  return new_sys
end

