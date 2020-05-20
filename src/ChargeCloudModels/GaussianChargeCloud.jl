function chargecloudsize(
        electric_field::Interpolations.Extrapolation{<:StaticVector{3}, 3}, 
        velocity_field::Interpolations.Extrapolation{<:StaticVector{3}, 3}, 
        cdm::AbstractChargeDriftModel{T}, 
        drift_path::Vector{CartesianPoint{T}}, 
        timestamps::Vector{T}, 
        charge::RealQuantity,
        energy::RealQuantity, 
        material;
        initial_size::RealQuantity = 0.1u"mm", 
        max_growth::RealQuantity = 0.05u"mm",
        temperature::RealQuantity = 77u"K"
    )::Vector{T} where {T <: SSDFloat}

    max_Δσ = to_internal_units.(internal_length_unit, max_growth)*10^9 #0.05 mm/ns
    max_ΔσΔt = to_internal_units.(internal_length_unit, 0.1u"mm")
    itemperature::T = ustrip(uconvert(u"K", temperature))

    n = length(timestamps)
    σ::Vector{T} = zeros(T, n)

    σ[1] = to_internal_units.(internal_length_unit, initial_size)
    previous_vel = norm(get_velocity_vector(velocity_field, drift_path[1]))

    ncharges =  to_internal_units(internal_energy_unit, energy) / to_internal_units(internal_energy_unit, material[:E_ionisation])
    repulsion_factor = ncharges * elementary_charge / (4*π*material.ϵ_r*ϵ0) * 0.68^3 #charge in 1σ, unitless
    
    for i in 2:n

        Δt = timestamps[i] - timestamps[i-1]
        σ[i] = σ[i-1]

        ####### acceleration part #######
        vel = norm(get_velocity_vector(velocity_field, drift_path[i]))
        σ[i] = σ[i-1] * vel / previous_vel

        ####### repulsion part #######
        Efield = get_velocity_vector(electric_field, drift_path[i])

        v_over_E = vel / norm(Efield)
        
        Efield_2 = Efield.*T(1.30) #arbitrary choice of 30% increase of E
        vel_2 = vel
        
        if charge < 0
            vel_2 = norm(getVe(Efield_2, cdm))
        else
            vel_2 = norm(getVh(Efield_2, cdm))
        end
        
        dv_dE = (vel_2 - vel) / (norm(Efield_2) - norm(Efield))
        
        #Δσ_rep = repulsion_factor * v_over_E / (σ[i]^2)
        Δσ_rep = repulsion_factor * dv_dE / (σ[i]^2)

        if(Δσ_rep > max_Δσ)
            Δσ_rep = max_Δσ
        end

        ####### diffusion part #######
        
        diff_coeff = kB * itemperature / elementary_charge * v_over_E # k*T/e is in siggen 0.67
        Δσ_dif = diff_coeff / σ[i]

        ####### final step #######

        Δσ_dt = Δσ_dif + Δσ_rep

        if(Δσ_dt > max_Δσ || (Δσ_dt * Δt) > max_ΔσΔt)
            σ[i] = σ[i] + sqrt(2 * diff_coeff * Δt + (σ[i]^2 * (3 * Δσ_rep * Δt))^(2/3))
            #println(i)
        else
            σ[i] = σ[i] + (Δσ_dif + Δσ_rep) * Δt
            #println(i, "\t", dv_dE, "\t", v_over_E)
        end
        
        previous_vel = vel

    end

    return σ

end
