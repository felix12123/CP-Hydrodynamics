
# structure for a hydrodynamic system
mutable struct HyDySys
  ρs::Vector{Float64}
  dx::Float64
  us::Vector{Float64}
  bound_cond::Symbol

  function HyDySys(ρs::Vector{Float64}, dx::Float64, us::Vector{Float64}=zeros(Float64, size(ρs)), bound_cond::Symbol=:periodic)
    valid_bound_conditions = [:periodic, :reflective] |> Tuple
    if !(bound_cond in valid_bound_conditions)
      error("boundary condition $bound_cond not valid. (valid conditions: $valid_bound_conditions)")
      return
    end
    new(ρs, dx, us, bound_cond)
  end
end



struct IntegratorScheme
  f::Function
  space_order::Int
  time_order::Int

  function IntegratorScheme(;f::Function=identity, space_order::Int=1, time_order::Int=1)
    new(f, space_order, time_order)
  end
end
