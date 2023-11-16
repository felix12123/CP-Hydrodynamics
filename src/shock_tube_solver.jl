using Statistics

function solve_shock_tube(sys::HyDySys, σ::Float64, a::Float64, t_ende::Float64)
    
    # Identifiziere wichtige Größen
    t       = 0
    l       = size(sys.us, 1)
    Δx      = sys.dx
    Δt      = σ * Δx / a
    ρ       = sys.ρs
    u       = sys.us
    ϵ       = sys.ϵs
    γ       = sys.γ
    x_order = 3


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
        if x > 0
            return 2 * x / (ρs[j+1]-ρs[j-1])
        else
            return 0
        end
        error("Fehler in Funktion: Δρ")
    end

    # (2.15)
    function ρ_adv(ρs, us, j, Δt, Δx)
        if us[j] > 0
          return ρs[j-1] + 0.5*(1-us[j]*Δt/Δx) * Δρ(ρs, j-1)
        else
          return ρs[j]   - 0.5*(1+us[j]*Δt/Δx) * Δρ(ρs, j)
        end
        error("Fehler in Funktion: ρ_adv")
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
        if x > 0
            return 2 * x / (us[j+1]-us[j-1])
        else
            return 0
        end
        error("Fehler in Funktion: Δu")
    end

    # (2.21)
    function mean_u(us, j)
        return mean([us[j], us[j+1]])
    end

    # (2.20)
    function u_adv(us, j, Δt, Δx)
        u_mean = mean_u(us, j) 
        if u_mean > 0 
            return us[j]   + 0.5 * (1 - u_mean * Δt/Δx) * Δu(us, j)
        else
            return us[j+1] - 0.5 * (1 + u_mean * Δt/Δx) * Δu(us, j+1)
        end
        error("Fehler in Funktion: u_adv")
    end

    # (2.19)
    function Fl(ρs, us, j, Δt, Δx)
        return mean([Fm(ρs, us, j, Δt, Δx), Fm(ρs, us, j+1, Δt, Δx)]) * u_adv(us, j, Δt, Δx)
    end

    # (2.23)
    function zwischenwert_u(ρs, us, j, Δt, Δx)
        return ((us[j]*mean([ρs[j], ρs[j-1]])) - Δt/Δx * (Fl(ρs, us, j, Δt, Δx) - Fl(ρs, us, j-1, Δt, Δx))) / mean([ρ_new(ρs, us, j-1, Δt, Δx), ρ_new(ρs, us, j, Δt, Δx)])
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
        error("Fehler in Funktion: Δϵ")
    end

    # (2.28)
    function ϵ_adv(ϵs, us, j, Δt, Δx) 
        if us[j] > 0
            return ϵs[j-1] + 0.5 * (1 - us[j] * Δt/Δx) * Δϵ(ϵs, j-1)
        else
            return ϵs[j]   - 0.5 * (1 + us[j] * Δt/Δx) * Δϵ(ϵs, j)
        end
        error("Fehler in Funktion: ϵ_adv")
    end

    # (2.27)
    function Fe(ρs, us, ϵs, j, Δt, Δx)
        return Fm(ρs, us, j, Δt, Δx) * ϵ_adv(ϵs, us, j, Δt, Δx)
    end

    # (2.30)
    function zwischenwert_ϵ(ρs, us, ϵs, j, Δt, Δx)
        return (ϵs[j]*ρs[j] + Δt/Δx * (Fe(ρs, us, ϵs, j+1, Δt, Δx) - Fe(ρs, us, ϵs, j, Δt, Δx))) / mean([ρ_new(ρs, us, j-1, Δt, Δx), ρ_new(ρs, us, j, Δt, Δx)])
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

    # (2.34)
    function u_new(ρs, us, ϵs, j, Δt, Δx, func_γ)
        return zwischenwert_u(ρs, us, j, Δt, Δx) - Δt * (1 / mean([ρ_new(ρs, us, j-1, Δt, Δx), ρ_new(ρs, us, j, Δt, Δx)])) * (zwischenwert_p(ρs, ϵs, func_γ, j) - zwischenwert_p(ρs, ϵs, func_γ, j - 1))/Δx 
    end
    # Nicht sicher, ob man den zwischenwert des Impulses aus der vorigen aufgabe nimmt oder nicht doch Gl. (2.36)

    # *******************
    # ***** Energie *****
    # *******************

    function ϵ_new(ρs, us, ϵs, j, Δt, Δx, func_γ)
        return zwischenwert_ϵ(ρs, us, ϵs, j, Δt, Δx) - Δt * (zwischenwert_p(ρs, ϵs, func_γ, j)/ρ_new(ρs, us, j, Δt, Δx)) * ((zwischenwert_u(ρs, us, j+1, Δt, Δx) - zwischenwert_u(ρs, us, j, Δt, Δx))/Δx)
    end
    # Nicht sicher, ob man den zwischenwert des Impulses aus der vorigen aufgabe nimmt oder nicht doch Gl. (2.36)

     
    
    
    # Beginne nun mit der eigentlichen Berechnung
    while t <= t_ende
        # Fertige Kopien der Systemgrößen an, damit bei der Berechnung nichts durcheinander kommt
        ρs = copy(ρ)
        us = copy(u)
        ϵs = copy(ϵ)

        
        # Füge Geister-Zellen basierend auf den Randbedingungen (2.39) ein
        if sys.bound_cond == :reflective
            us = vcat([-us[2], 0],     us, [0,      -us[end-1]])
            ρs = vcat([ρs[2],  ρs[1]], ρs, [ρs[end], ρs[end-1]])
            ϵs = vcat([ϵs[2],  ϵs[1]], ϵs, [ϵs[end], ϵs[end-1]])
        elseif sys.bound_cond == :periodic
            ρs = vcat(ρs[end-x_order+1:end], ρs, ρs[1:x_order])
            us = vcat(us[end-x_order+1:end], us, us[1:x_order])
            ϵs = vcat(ϵs[end-x_order+1:end], ϵs, ϵs[1:x_order])
        else
           error("Boundary condition not valid: $(sys.bound_cond)")
           return false
        end
        
        # Fertige zusätzlich nochmal Kopien an, die nicht geupdated werden, damit
        # der zweite Zwischen- schritt der Berechnung korrekt durchgeführt werden kann.
        # So wird die Reihenhfolge bei der Berechnung des Advektionsschritts irrelevant.
        ρs_not_updated = copy(ρs)
        ϵs_not_updated = copy(ϵs)

        # Berechnung des Advektionsschritts, die updates erfolgen nach (2.31)
        for j in (x_order+1):(l + x_order)
            ρs[j] = ρ_new(ρs, us, j, Δt, Δx)
            us[j] = zwischenwert_u(ρs_not_updated, us,     j, Δt, Δx)
            ϵs[j] = zwischenwert_ϵ(ρs_not_updated, us, ϵs, j, Δt, Δx)
        end
        
        # Füge Geister-Zellen basierend auf den Randbedingungen (2.39) ein
        if sys.bound_cond == :reflective
            us = vcat([-us[2], 0],     us[3:end-2], [0,      -us[end-1]])
            ρs = vcat([ρs[2],  ρs[1]], ρs[3:end-2], [ρs[end], ρs[end-1]])
            ϵs = vcat([ϵs[2],  ϵs[1]], ϵs[3:end-2], [ϵs[end], ϵs[end-1]])
        elseif sys.bound_cond == :periodic
            ρs = vcat(ρs[end-x_order+1:end], ρs[(x_order+1):(end-x_order)], ρs[1:x_order])
            us = vcat(us[end-x_order+1:end], us[(x_order+1):(end-x_order)], us[1:x_order])
            ϵs = vcat(ϵs[end-x_order+1:end], ϵs[(x_order+1):(end-x_order)],ϵs[1:x_order])
        else
           error("Boundary condition not valid: $(sys.bound_cond)")
           return false
        end

        # Berechnung der Kräfte und Druckarbeit, hier werden u und ϵ geupdated
        for j in (x_order+1):(l + x_order)
            us[j] = u_new(ρs, us, ϵs, j, Δt, Δx, γ)
            ϵs[j] = ϵ_new(ρs, us, ϵs, j, Δt, Δx, γ)
        end

        us = us[(x_order+1):(end-x_order)]
        ρs = ρs[(x_order+1):(end-x_order)]
        ϵs = ϵs[(x_order+1):(end-x_order)]

        # Update Zeitschritt
        t += Δt
    end

    return HyDySys(ρ, Δx, u, :reflective, ϵ, γ)
end