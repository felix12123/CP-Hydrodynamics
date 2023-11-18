# σ = a dt/dx   =>   dt = dx σ / a
function solve_lin_adv(sys::HyDySys, σ::Float64, a::Float64, t_end::Float64)
  space_order::Int = 2
  t::Float64       = 0.0
  N::Int           = size(sys.us, 1)
  dt::Float64      = σ * sys.dx / a
  dx::Float64      = sys.dx
  γ::Float64       = sys.γ
  
  ρs::Vector{Float64} = copy(sys.ρs)
  us::Vector{Float64} = copy(sys.us)
  ϵs::Vector{Float64} = copy(sys.ϵs)
  
  
  
  # Define Functions that are needed for better readability of the time step update
  function Δρ(ρs::Vector{Float64}, j::Int)::Float64
    x = (ρs[j+1] - ρs[j])*(ρs[j] - ρs[j-1]) 
    if x > 0
      2 * x / (ρs[j+1] - ρs[j-1])
    else
      0
    end
  end

  function ρ_adv(ρs::Vector{Float64}, us::Vector{Float64}, j::Int, dt_dx::Float64)::Float64
    # return ρs[j-1]
    if us[j] > 0
      return ρs[j-1] + 1/2*(1-us[j]*dt_dx) * Δρ(ρs, j-1)
    else
      return ρs[j]   - 1/2*(1+us[j]*dt_dx) * Δρ(ρs, j)
    end
  end
  
  function Fm(ρs::Vector{Float64}, us::Vector{Float64}, j::Int, dt_dx::Float64)::Float64
    ρ_adv(ρs, us, j, dt_dx) * us[j]
  end
  counter = 0
  # Start the actual simulation
  while t < t_end
    counter += 1
    # add ghost cells, according to boundary condition
    if sys.bound_cond == :periodic
      ρs = vcat(ρs[end-space_order+1:end], ρs, ρs[1:space_order])
      us = vcat(us[end-space_order+1:end], us, us[1:space_order])
    else # might incorporate other boundary conditions later
      error("Boundary condition not valid: $(sys.bound_cond)")
      return false
    end

    # copy current state, to avoid aliasing
    ρs_copy = copy(ρs)
    us_copy = copy(us)

    # update every cell
    for j in space_order+1:N+space_order
      ρs[j] = ρs_copy[j] - dt/dx * (Fm(ρs_copy, us_copy, j+1, dt/dx) - Fm(ρs_copy, us_copy, j, dt/dx))
    end
    t += dt

    # remove ghost cells. they either have to be renewed in the next step,
    # or they have to be cut before the system is returned
    ρs = ρs[space_order+1:end-space_order]
    us = us[space_order+1:end-space_order]
  end
  println(ρs)
  new_sys = HyDySys(ρs, dx, us, sys.bound_cond, ϵs, γ)
  return new_sys
end