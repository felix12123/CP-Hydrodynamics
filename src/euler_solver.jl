function solve_euler(N, γ)

    # *****************************
	# ***** Advektionsschritt *****
	# *****************************

	# ******************
	# ***** Dichte *****
	# ******************

    function ρ_advection!(Δt, Δx, ρ, u, ϵ, mass_flux, Δρ, ρ_adv, N)

        # (2.16) - Berechnung der Δρ
        for i in 3:(N-2)
            cond = (ρ[i+1]-ρ[i])*(ρ[i]-ρ[i-1])
            if cond > 0
                Δρ[i] = 2*cond/(ρ[i+1] - ρ[i-1])
            else
                Δρ[i] = 0
            end
        end

        # (2.15) - Berechnung ρ_adv
        for i in 3:(N-2)
            if u[i] > 0
                ρ_adv[i] = ρ[i-1] + 0.5*(1-u[i]*Δt/Δx)*Δρ[i-1]
            else
                ρ_adv[i] = ρ[i] + 0.5*(1+u[i]*Δt/Δx)*Δρ[i]
            end
        end

        # (2.14) - Berechnung Massenfluss
        for i in 3:(N-2)
            mass_flux[i] = ρ_adv[i]*u[i]
        end

        # (2.13)
        for i in 3:(N-2)
            ρ[i] = ρ[i] -  Δt/ Δx*(mass_flux[i+1] - mass_flux[i])
        end

    end



    # **********************************
	# ***** Impuls/Geschwindigkeit *****
	# **********************************

    function u_advection!(Δt, Δx, ρ, ρ0, u, ϵ, mass_flux, momentum_flux, Δu, u_adv, N)

        # (2.22) - Berechne Δu
        for i in 3:(N-2)
            cond = (u[i+1]-u[i])*(u[i]-u[i-1])
            if cond > 0
                Δu[i] = 2*cond/(u[i+1]-u[i-1])
            else
                Δu[i] = 0
            end
        end

        # (2.20) - Berechne u_adv
        for i in 3:(N-2)
            ubar = 0.5*(u[i]+u[i+1])
            if ubar > 0
                u_adv[i] = u[i] + 0.5*(1-ubar*Δt/Δx)*Δu[i]
            else
                u_adv[i] = u[i+1] - 0.5*(1+ubar*Δt/Δx)*Δu[i+1] 
            end
        end

        # (2.19) - Berechne Impulsfluss
        for i in 3:(N-2)
            momentum_flux[i] = 0.5*(mass_flux[i]+mass_flux[i+1])*u_adv[i]
        end

        # (2.23) - Berechne Zwischenwert von u
        for i in 3:(N-2)
            ρbar0 = 0.5*(ρ0[i-1] + ρ0[i])
            ρbar  = 0.5*(ρ[i-1] + ρ[i])
            u[i] = 1/ρbar*(u[i]*ρbar0-Δt/Δx*(momentum_flux[i]-momentum_flux[i-1]))
        end

    end



    # *******************
	# ***** Energie *****
	# *******************

    function ϵ_advection!(Δt, Δx, ρ, ρ0, u, ϵ, mass_flux, energy_flux, Δϵ, ϵ_adv, N)
    
        # (2.29) - Berechne Δϵ
        for i in 3:(N-2)
            cond = (ϵ[i+1]-ϵ[i])*(ϵ[i]-ϵ[i-1])
            if cond > 0
                Δϵ[i] = 2*cond/(ϵ[i+1]-ϵ[i-1])
            else
                Δϵ[i] = 0
            end
        end

        # (2.28) - Berechne ϵ_adv
        for i in 3:(N-2)
            if u[i] > 0
                ϵ_adv[i] = ϵ[i-1] + 0.5*(1-u[i]*Δt/Δx)*Δϵ[i-1]
            else
                ϵ_adv[i] = ϵ[i]   - 0.5*(1+u[i]*Δt/Δx)*Δϵ[i]
            end
        end

        # (2.27) - Berechne den Energiefluss
        for i in 3:(N-2)
            energy_flux[i] = mass_flux[i] * ϵ_adv[i]
        end

        # (2.30) - Berechne Zwischenwert von ϵ
        for i in 3:(N-2)
            ϵ[i] = 1/ρ[i] * (ϵ[i]*ρ0[i] - Δt/Δx *(energy_flux[i+1] - energy_flux[i]))
        end

    end



    function calculate_pressure!(p, ρ, ϵ, N, γ)

        # (2.36) - Berechne Druck
        for i in 3:(N-2)
            p[i] = (γ-1)*ρ[i]*ϵ[i]
        end

    end



    function calculate_Temperature!(T, ϵ, N, γ)

        # Berechne Temperatur
        for i in 3:(N-2)
            T[i] = ϵ[i]*(γ-1)
        end

    end



    function update_u!(Δt, Δx, ρ, u, p, N)

        # (2.34) - Berechne Endwert von u
        for i in 3:(N-2)
            ρbar = 0.5*(ρ[i-1] + ρ[i])
            u[i] = u[i] -Δt/Δx*(p[i] - p[i-1])/ρbar
        end

    end



    function update_ϵ!(Δt, Δx, ρ, u_old, ϵ, p, N)

        # (2.38) - Berechne Endwert von ϵ
        for i in 3:(N-2)
            ϵ[i] = ϵ[i] -Δt/Δx*p[i]/ρ[i] * (u_old[i+1] - u_old[i])
        end

    end



	# *****************************
	# ***** Starte Berechnung *****
	# *****************************

    # Parameter aus der Aufgabenstellunge
    x   = range(0,1,N)
    Δx  = x[2]-x[1]
    x_0 = 0.5
    Δt  = 0.001
    t_current = 0.0
    t_end     = 0.228

    # Initialisiere alle Arrays
    ρ             = zeros(Float64, N)
    ρ0            = zeros(Float64, N)
    ϵ             = zeros(Float64, N)
    pressure      = zeros(Float64, N)
    temperature   = zeros(Float64, N)
    u             = zeros(Float64, N+1)
    u_old         = zeros(Float64, N+1)
    mass_flux     = zeros(Float64, N+1)
    energy_flux   = zeros(Float64, N+1)
    momentum_flux = zeros(Float64, N+1)
    Δu            = zeros(Float64, N+1)
    u_adv         = zeros(Float64, N+1)
    Δϵ            = zeros(Float64, N)
    ϵ_adv         = zeros(Float64, N)
    Δρ            = zeros(Float64, N)
    ρ_adv         = zeros(Float64, N)

    # Generiere Startsystem
    for i in 1:N
		if x[i] <= x_0
			pressure[i] = 1.0
			ρ[i]        = 1.0
			ϵ[i]        = 2.5
		else
			pressure[i] = 0.1
			ρ[i]        = 0.125
			ϵ[i]        = 2.0
		end
	end

    # Beginne mit Algorithmus
    while t_current < t_end

        ρ0 = deepcopy(ρ)

        calculate_pressure!(pressure, ρ, ϵ, N, γ)

        ρ_advection!(Δt, Δx, ρ, u, ϵ, mass_flux, Δρ, ρ_adv, N)

        u_advection!(Δt, Δx, ρ, ρ0, u, ϵ, mass_flux, momentum_flux, Δu, u_adv, N)
            
        ϵ_advection!(Δt, Δx, ρ, ρ0, u, ϵ, mass_flux, energy_flux, Δϵ, ϵ_adv, N)

        calculate_pressure!(pressure, ρ, ϵ, N, γ)

        u_old = deepcopy(u)

        update_u!(Δt, Δx, ρ, u, pressure, N)

        update_ϵ!(Δt, Δx, ρ, u_old, ϵ, pressure, N)

        calculate_Temperature!(temperature, ϵ, N, γ)

        t_current += Δt
    end

    plot(x, u[1:end-1], title="", xlabel=L"x", ylabel=L"u", dpi=300, linewidth=1, linealpha=0.2, linecolor = :black, label="")
    scatter!(x, u[1:end-1], marker=:xcross, markersize=2, markerstrokewidth=2, markercoloer= :orange,label="")
    savefig("media/A2_u_228")
        
    plot(x, ρ, title="", xlabel=L"x", ylabel=L"\rho", dpi=300, linewidth=1, linealpha=0.2, linecolor = :black, label="")
    scatter!(x, ρ, marker=:xcross, markersize=2, markerstrokewidth=2, markercoloer= :orange,label="")
    savefig("media/A2_rho_228")
        
    plot(x, temperature, title="", xlabel=L"x", ylabel=L"T", dpi=300, linewidth=1, linealpha=0.2, linecolor = :black, label="")
    scatter!(x, temperature, marker=:xcross, markersize=2, markerstrokewidth=2, markercoloer= :orange,label="")
    savefig("media/A2_T_228")

    plot(x, ϵ, title="", xlabel=L"x", ylabel=L"\epsilon", dpi=300, linewidth=1, linealpha=0.2, linecolor = :black, label="")
    scatter!(x, ϵ, marker=:xcross, markersize=2, markerstrokewidth=2, markercoloer= :orange,label="")
    savefig("media/A2_epsilon_228")

end