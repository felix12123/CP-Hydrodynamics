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