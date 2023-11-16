# structure for a hydrodynamic system
mutable struct HyDySys
  ρs::Vector{Float64}
  dx::Float64
  us::Vector{Float64}
  bound_cond::Symbol
  ϵs::Vector{Float64}
  γ::Float64

  function HyDySys(ρs::Vector{Float64}, dx::Float64, us::Vector{Float64}=zeros(Float64, size(ρs)), bound_cond::Symbol=:periodic, ϵs::Vector{Float64}=ones(Float64, size(ρs)), γ=0.5)
    valid_bound_conditions = [:periodic, :reflective] |> Tuple
    if !(bound_cond in valid_bound_conditions)
      error("boundary condition $bound_cond not valid. (valid conditions: $valid_bound_conditions)")
      return
    end
    new(ρs, dx, us, bound_cond, ϵs, γ)
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
