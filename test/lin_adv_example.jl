function small_test()
	T  = 0.0
	N  = 2<<8
	dx = 1.0
	a  = 1.0
	σ=0.8

	pack_size = ceil(Int, N/4)
	ρs = vcat(zeros(Float64, floor(Int, (N-pack_size)/2)), ones(Float64, pack_size), zeros(Float64, ceil(Int, (N-pack_size)/2)))
	us = ones(Float64, N) .* 2.0

	sys = HyDySys(ρs, dx, us, :periodic)
	new_sys = solve_lin_adv(sys, σ, a, T)

	# display([sys, new_sys])
	display(new_sys)
end