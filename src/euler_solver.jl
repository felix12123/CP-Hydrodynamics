using Statistics

function solve_euler(sys::HyDySys, σ::Float64, a::Float64, t_ende::Float64)
    
    # Identifiziere wichtige Größen
    t       = 0.0
    l       = size(sys.us, 1)
    Δx      = sys.dx
    Δt      = σ * Δx / a
    ρs      = copy(sys.ρs)
    us      = copy(sys.us)
    ϵs      = copy(sys.ϵs)
    γ       = sys.γ
    Ordnung = 2

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
        critical = ρs[j+1]-ρs[j-1]
        if x > 0
            if critical == 0
                return 2 * x / 0.000001
            elseif critical != 0
                return 2 * x / (ρs[j+1]-ρs[j-1])
            end
        else
            return 0
        end
    end

    # (2.15)
    function ρ_adv(ρs, us, j, Δt, Δx)
        if us[j] > 0
            return ρs[j-1] + (1/2)*(1-us[j]*Δt/Δx) * Δρ(ρs, j-1)
        else
            return ρs[j]   - (1/2)*(1+us[j]*Δt/Δx) * Δρ(ρs, j)
        end
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
        critical = us[j+1]-us[j-1]
        if x > 0
            if critical == 0
                return 2 * x / 0.000001
            elseif critical != 0
                return 2 * x / critical
            end
        else
            return 0
        end
        error("Fehler in Funktion: Δu")
    end

    # (2.20)
    function u_adv(us, j, Δt, Δx)
        u_mean = mean([us[j], us[j+1]]) 
        if u_mean > 0 
            return us[j]   + (1/2) * (1 - u_mean * Δt/Δx) * Δu(us, j)
        else
            return us[j+1] - (1/2) * (1 + u_mean * Δt/Δx) * Δu(us, j+1)
        end
        error("Fehler in Funktion: u_adv")
    end

    # (2.19)
    function Fl(ρs, us, j, Δt, Δx)
        return mean([Fm(ρs, us, j, Δt, Δx), Fm(ρs, us, j+1, Δt, Δx)]) * u_adv(us, j, Δt, Δx)
    end

    # (2.23)
    function zwischenwert_u(ρs, ρs_new, us, j, Δt, Δx)
        critical = mean([ρs_new[j-1], ρs_new[j]])
        if critical == 0
            return ((us[j]*mean([ρs[j], ρs[j-1]])) - Δt/Δx * (Fl(ρs, us, j, Δt, Δx) - Fl(ρs, us, j-1, Δt, Δx))) / 0.000001
        else
            return ((us[j]*mean([ρs[j], ρs[j-1]])) - Δt/Δx * (Fl(ρs, us, j, Δt, Δx) - Fl(ρs, us, j-1, Δt, Δx))) / critical
        end
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
            return ϵs[j-1] + (1/2) * (1 - us[j] * Δt/Δx) * Δϵ(ϵs, j-1)
        else
            return ϵs[j]   - (1/2) * (1 + us[j] * Δt/Δx) * Δϵ(ϵs, j)
        end
        error("Fehler in Funktion: ϵ_adv")
    end

    # (2.27)
    function Fe(ρs, us, ϵs, j, Δt, Δx)
        return Fm(ρs, us, j, Δt, Δx) * ϵ_adv(ϵs, us, j, Δt, Δx)
    end

    # (2.30)
    function zwischenwert_ϵ(ρs, ρs_new, us, ϵs, j, Δt, Δx)
        critical = mean([ρs_new[j-1], ρs_new[j]])
        if critical == 0
            return (ϵs[j]*ρs[j] + Δt/Δx * (Fe(ρs, us, ϵs, j+1, Δt, Δx) - Fe(ρs, us, ϵs, j, Δt, Δx))) / 0.000001
        else
            return (ϵs[j]*ρs[j] + Δt/Δx * (Fe(ρs, us, ϵs, j+1, Δt, Δx) - Fe(ρs, us, ϵs, j, Δt, Δx))) / critical
        end
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

    function u_new(ρs_new, us_zwischen, ϵs_zwischen, j ,Δt, Δx, func_γ)
        critical = mean([ρs_new[j-1], ρs_new[j]])
        if critical == 0
            return us_zwischen[j] - Δt * (1/0.00001) * ((func_γ-1)*(ρs_new[j]ϵs_zwischen[j]-ρs_new[j-1]ϵs_zwischen[j-1])/Δx)
        else
            return us_zwischen[j] - Δt * (1/critical) * ((func_γ-1)*(ρs_new[j]ϵs_zwischen[j]-ρs_new[j-1]ϵs_zwischen[j-1])/Δx)
        end
    end

    # *******************
    # ***** Energie *****
    # *******************   
    
    function ϵ_new(ρs_new, us_zwischen, ϵs_zwischen, j ,Δt, Δx, func_γ)
        critical = ρs_new[j+1]
        if critical == 0
            return ϵs_zwischen[j] - Δt * (((func_γ-1)*ρs_new[j]*ϵs_zwischen[j])/0.000001) * ((us_zwischen[j+1]-us_zwischen[j])/Δx)
        elseif critical != 0
            return ϵs_zwischen[j] - Δt * (((func_γ-1)*ρs_new[j]*ϵs_zwischen[j])/ρs_new[j+1]) * ((us_zwischen[j+1]-us_zwischen[j])/Δx)
        end
    end

    # Beginne nun mit der eigentlichen Berechnung
    while t < t_ende
        # Füge Geister-Zellen basierend auf den Randbedingungen (2.39) ein
        if sys.bound_cond == :reflective
            us = vcat([-us[2], 0],     us, [0,      -us[end-1]])
            ρs = vcat([ ρs[2], ρs[1]], ρs, [ρs[end], ρs[end-1]])
            ϵs = vcat([ ϵs[2], ϵs[1]], ϵs, [ϵs[end], ϵs[end-1]])
        elseif sys.bound_cond == :periodic
            ρs = vcat(ρs[end-Ordnung+1:end], ρs, ρs[1:Ordnung])
            us = vcat(us[end-Ordnung+1:end], us, us[1:Ordnung])
            ϵs = vcat(ϵs[end-Ordnung+1:end], ϵs, ϵs[1:Ordnung])
        else
            error("Boundary condition not valid: $(sys.bound_cond)")
            return false
        end
        # Fertige Kopien der Systemgrößen an, damit bei der Berechnung nichts durcheinander kommt
        ρs_copy = copy(ρs)
        us_copy = copy(us)
        ϵs_copy = copy(ϵs)
        
        # Fertige zusätzlich nochmal Kopien an, die nicht geupdated werden, damit
        # der zweite Zwischen- schritt der Berechnung korrekt durchgeführt werden kann.
        # So wird die Reihenhfolge bei der Berechnung des Advektionsschritts irrelevant.
        ρs_not_updated = copy(ρs_copy)
        ϵs_not_updated = copy(ϵs_copy)

        # Berechnung des Advektionsschritts, die updates erfolgen nach (2.31):
        # ρ → ρ^{n+1}
        # u → u^1
        # ϵ → ϵ^1
        # Berechne zuerst die neue Dichte des Systems:
        for j in (Ordnung+1):(l + Ordnung - 1)
            ρs[j] = ρ_new(ρs_copy, us_copy, j, Δt, Δx)

        end

        # Berechne nun die Zwischenwerte der Geschwindigkeit und Energie
        for j in (Ordnung+1):(l + Ordnung - 1)
            # If-Schleife damit [j=3...N+1] eingehalten wird
            if sys.bound_cond == :reflective
                if (j >= (Ordnung+2))
                    us[j] = zwischenwert_u(ρs_not_updated, ρs, us_copy, j, Δt, Δx)
                end
                ϵs[j] = zwischenwert_ϵ(ρs_not_updated, ρs, us_copy, ϵs_not_updated, j, Δt, Δx) 
            end
        end

        #println("Dichte zu t=", t, ": " , ρs)
        #println("Zwischen-Geschwindigkeit:", us)
        #println("Zwischen-Energie:", ϵs)

        # Füge Geister-Zellen basierend auf den Randbedingungen (2.39) ein
        if sys.bound_cond == :reflective
            us = vcat([-us[2], 0],     us[3:end-2], [0,      -us[end-1]])
            ρs = vcat([ ρs[2], ρs[1]], ρs[3:end-2], [ρs[end], ρs[end-1]])
            ϵs = vcat([ ϵs[2], ϵs[1]], ϵs[3:end-2], [ϵs[end], ϵs[end-1]])
        elseif sys.bound_cond == :periodic
            nothing
        else
            error("Boundary condition not valid: $(sys.bound_cond)")
            return false
        end

        # Berechnung der Kräfte und Druckarbeit, hier werden u und ϵ geupdated
        # u → u^1 → u^{n+1}
        # ϵ → ϵ^1 → ϵ^{n+1}
        if sys.bound_cond == :reflective    
            for j in (Ordnung + 1):(l + Ordnung - 1)
                # If-Schleife damit [j=3...N+1] eingehalten wird
                if j >= (Ordnung+2)
                    us[j] = u_new(ρs, us, ϵs, j, Δt, Δx, γ)
                end
                ϵs[j] = ϵ_new(ρs, us, ϵs, j, Δt, Δx, γ)
            end
        end

        #println("Update-Geschwindigkeit:", us)
        #println("Update-Energie:", ϵs)

        # Entferne Giisterzellen
        us = us[(Ordnung+1):(end-Ordnung)]
        ρs = ρs[(Ordnung+1):(end-Ordnung)]
        ϵs = ϵs[(Ordnung+1):(end-Ordnung)]
        

        # Update Zeitschritt
        t += Δt
    end

    return HyDySys(ρs, Δx, us, sys.bound_cond, ϵs, γ)
end