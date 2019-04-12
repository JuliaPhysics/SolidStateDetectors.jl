abstract type AbstractSSDSetup{T <: SSDFloat} end

mutable struct SSDSetup{T <: SSDFloat} <: AbstractSSDSetup{T}
    detector::SolidStateDetector{T}
    ρ::ChargeDensity{T}
    ϵ::DielectricDistribution{T}
    point_types::PointTypes{T}
    electric_potential::ElectricPotential{T}  
    weighting_potentials::Dict{Int, WeightingPotential{T}}
    electric_field::ElectricField{T}

    charge_drift_model::AbstractChargeDriftModel{T}

    interpolated_electron_drift_field
    interpolated_hole_drift_field

    SSDSetup{T}() where {T <: SSDFloat} = new{T}()
end

function SSDSetup(detector::SolidStateDetector{T})::SSDSetup{T} where {T <: SSDFloat}
    setup::SSDSetup{T} = SSDSetup{T}()
    setup.detector = detector
    setup.weighting_potentials = Dict{Int, WeightingPotential{T}}()
    return setup
end

function SSDSetup{T}(config_file::AbstractString)::SSDSetup{T} where{T <: SSDFloat}
    return SSDSetup( SolidStateDetector{T}(config_file) )
end

function apply_initial_state!(setup::SSDSetup{T})::Nothing where {T <: SSDFloat}
    init_setup = SolidStateDetectors.PotentialSimulationSetup(setup.detector);

    setup.ρ = ChargeDensity(init_setup.ρ, init_setup.grid)
    setup.ϵ = DielectricDistribution(init_setup.ϵ, init_setup.grid)
    setup.point_types = PointTypes(init_setup.pointtypes, init_setup.grid)
    setup.electric_potential = ElectricPotential(init_setup.potential, init_setup.grid)

    nothing
end

function calculate_electric_potential!(setup::SSDSetup{T}, args...; kwargs...)::Nothing where {T <: SSDFloat}
    ep = calculate_electric_potential(setup.detector, args...; kwargs...)

    setup.ρ = ChargeDensity(ep.ρ, ep.grid)
    setup.ϵ = DielectricDistribution(ep.ϵ, ep.grid)
    setup.point_types = PointTypes(ep.pointtypes, ep.grid)
    setup.electric_potential = ElectricPotential(ep.potential, ep.grid)

    nothing
end

function calculate_weighting_potential!(setup::SSDSetup{T}, contact_id::Int, args...; n_points_in_φ::Union{Missing, Int} = missing, kwargs...)::Nothing where {T <: SSDFloat}
    S = SSD.get_coordinate_system(setup.detector)
    if S == :Cylindrical && setup.detector.cyclic == T(0)
        if ismissing(n_points_in_φ)
            @info "\tIn weighing potential calculation: Keyword `n_points_in_φ` not set.\n\t\tDefault is `n_points_in_φ = 18`. 2D field will be extended to 18 points in φ."
            n_points_in_φ = 18
        else
            if !(n_points_in_φ > 1 && iseven(n_points_in_φ))
                @info "\tIn weighing potential calculation: Keyword `n_points_in_φ` is $(n_points_in_φ) but must be even and larger than 1.\n\t\t`n_points_in_φ` is now set to 18. 2D field will be extended to 18 points in φ."
                n_points_in_φ = 18
            end
        end
    end
    wps = calculate_weighting_potential(setup.detector, contact_id, args...; kwargs...)
    if S == :Cylindrical && size(wps.potential, 2) == 1 && !ismissing(n_points_in_φ)
        wp = WeightingPotential(wps, n_points_in_φ = n_points_in_φ)
    else
        wp = WeightingPotential(wps)
    end
    push!(setup.weighting_potentials, (contact_id => wp))
    nothing
end


function calculate_electric_field!(setup::SSDSetup{T}, args...; n_points_in_φ::Union{Missing, Int} = missing, kwargs...)::Nothing where {T <: SSDFloat}
    S = SSD.get_coordinate_system(setup.detector)
    e_pot, point_types = if S == :Cylindrical && setup.detector.cyclic == T(0) # 2D, only one point in φ
        if ismissing(n_points_in_φ)
            @info "\tIn electric field calculation: Keyword `n_points_in_φ` not set.\n\t\tDefault is `n_points_in_φ = 18`. 2D field will be extended to 18 points in φ."
            n_points_in_φ = 18
        else
            if !(n_points_in_φ > 1 && iseven(n_points_in_φ))
                @info "\tIn electric field calculation: Keyword `n_points_in_φ` is $(n_points_in_φ) but must be even and larger than 1.\n\t\t`n_points_in_φ` is now set to 18. 2D field will be extended to 18 points in φ."
                n_points_in_φ = 18
            end
        end
        get_2π_potential(setup.electric_potential, n_points_in_φ = n_points_in_φ),
        get_2π_potential(setup.point_types,  n_points_in_φ = n_points_in_φ);
    else
        setup.electric_potential,
        setup.point_types
    end
    setup.electric_field = get_electric_field_from_potential(e_pot, point_types);
    nothing
end

function set_charge_drift_model!(setup::SSDSetup{T}, charge_drift_model::AbstractChargeDriftModel)::Nothing where {T <: SSDFloat}
    setup.charge_drift_model = charge_drift_model
    nothing
end

function apply_charge_drift_model!(setup::SSDSetup{T})::Nothing where {T <: SSDFloat}
    setup.interpolated_electron_drift_field = get_interpolated_drift_field(
        get_electron_drift_field(setup.electric_field.data, setup.charge_drift_model), setup.electric_field.grid
    )
    setup.interpolated_hole_drift_field = get_interpolated_drift_field(
        get_hole_drift_field(setup.electric_field.data, setup.charge_drift_model), setup.electric_field.grid
    )
    nothing
end

@inline function drift_charges( setup::SSDSetup{T}, starting_positions::Vector{CartesianPoint{T}}, 
                                Δt::T = T(1f-9), n_steps::Int = 2000 )::DriftPath{T} where {T <: SSDFloat}
    return drift_charges(   setup.detector, setup.electric_potential.grid, starting_positions, 
                            setup.interpolated_electron_drift_field, setup.interpolated_hole_drift_field,
                            Δt = Δt, n_steps = n_steps)::DriftPath{T}
end

