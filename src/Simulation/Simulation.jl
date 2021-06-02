abstract type AbstractSimulation{T <: SSDFloat} end

"""
    mutable struct Simulation{T <: SSDFloat, CS <: AbstractCoordinateSystem} <: AbstractSimulation{T}

Collection of all parts of a Simulation of a Solid State Detector.
"""
mutable struct Simulation{T <: SSDFloat, CS <: AbstractCoordinateSystem} <: AbstractSimulation{T}
    config_dict::Dict
    input_units::NamedTuple
    medium::NamedTuple # this should become a struct at some point
    detector::Union{SolidStateDetector{T}, Missing}
    world::World{T, 3, CS}
    q_eff_imp::Union{EffectiveChargeDensity{T}, Missing} # Effective charge coming from the impurites of the semiconductors
    q_eff_fix::Union{EffectiveChargeDensity{T}, Missing} # Fixed charge coming from fixed space charges, e.g. charged up surface layers
    ϵ_r::Union{DielectricDistribution{T}, Missing}
    point_types::Union{PointTypes{T}, Missing}
    electric_potential::Union{ElectricPotential{T}, Missing}
    weighting_potentials::Vector{Any}
    electric_field::Union{ElectricField{T}, Missing}
    electron_drift_field::Union{ElectricField{T}, Missing}
    hole_drift_field::Union{ElectricField{T}, Missing}
end

function Simulation{T,CS}() where {T <: SSDFloat, CS <: AbstractCoordinateSystem}
    Simulation{T, CS}(
        Dict(),
        default_unit_tuple(),
        material_properties[materials["vacuum"]],
        SolidStateDetector{T}(),
        World(CS,(T(0),T(1),T(0),T(1),T(0),T(1))),
        missing,
        missing,
        missing,
        missing,
        missing,
        [missing],
        missing,
        missing,
        missing
    )
end

get_precision_type(::Simulation{T}) where {T} = T
get_coordinate_system(::Simulation{T, CS}) where {T, CS} = CS

function NamedTuple(sim::Simulation{T}) where {T <: SSDFloat}
    wps_strings = AbstractString[]
    for contact in sim.detector.contacts
        if !ismissing(sim.weighting_potentials[contact.id])
            push!(wps_strings, "WeightingPotential_$(contact.id)")
        end
    end
    wps_syms = Symbol.(wps_strings)
    nt = (
        detector_json_string = NamedTuple(sim.detector.config_dict),
        electric_potential = NamedTuple(sim.electric_potential),
        q_eff_imp = NamedTuple(sim.q_eff_imp),
        q_eff_fix = NamedTuple(sim.q_eff_fix),
        ϵ_r = NamedTuple(sim.ϵ_r),
        point_types = NamedTuple(sim.point_types),
        electric_field = NamedTuple(sim.electric_field),
        electron_drift_field = NamedTuple(sim.electron_drift_field),
        hole_drift_field = NamedTuple(sim.hole_drift_field)
    )
    if length(wps_strings) > 0
        nt = merge(nt, (weighting_potentials = NamedTuple{Tuple(wps_syms)}(NamedTuple.( skipmissing(sim.weighting_potentials))),))
    end
    return nt
end
Base.convert(T::Type{NamedTuple}, x::Simulation) = T(x)

function Simulation(nt::NamedTuple)
    ep = ElectricPotential(nt.electric_potential)
    T = eltype(ep.data)
    det = SolidStateDetector{T}( Dict(nt.detector_json_string) )
    sim = Simulation( det )
    if !ismissing(nt.electric_potential) sim.electric_potential = ep end
    if !ismissing(nt.q_eff_imp) sim.q_eff_imp = EffectiveChargeDensity(nt.q_eff_imp) end
    if !ismissing(nt.q_eff_fix) sim.q_eff_fix = EffectiveChargeDensity(nt.q_eff_fix) end
    if !ismissing(nt.ϵ_r) sim.ϵ_r = DielectricDistribution(nt.ϵ_r) end
    if !ismissing(nt.point_types) sim.point_types = PointTypes(nt.point_types) end
    if !ismissing(nt.electric_field) sim.electric_field = ElectricField(nt.electric_field) end
    sim.weighting_potentials = [missing for contact in sim.detector.contacts]
    if !ismissing(nt.weighting_potentials)
        for contact in sim.detector.contacts
            if !ismissing(values(nt.weighting_potentials[contact.id])[1])
                sim.weighting_potentials[contact.id] = WeightingPotential(nt.weighting_potentials[contact.id])
            end
        end
    end
    if !ismissing(nt.electron_drift_field) sim.electron_drift_field = ElectricField(nt.electron_drift_field) end
    if !ismissing(nt.hole_drift_field) sim.hole_drift_field = ElectricField(nt.hole_drift_field) end
    return sim
end
Base.convert(T::Type{Simulation}, x::NamedTuple) = T(x)




function println(io::IO, sim::Simulation{T}) where {T <: SSDFloat}
    println(typeof(sim), " - Coordinate system: ", get_coordinate_system(sim))
    println("  Environment Material: $(sim.medium.name)")
    println("  Detector: $(sim.detector.name)")
    println("  Electric potential: ", !ismissing(sim.electric_potential) ? size(sim.electric_potential) : missing)
    println("  Charge density: ", !ismissing(sim.q_eff_imp) ? size(sim.q_eff_imp) : missing)
    println("  Fix Charge density: ", !ismissing(sim.q_eff_fix) ? size(sim.q_eff_fix) : missing)
    println("  Dielectric distribution: ", !ismissing(sim.ϵ_r) ? size(sim.ϵ_r) : missing)
    println("  Point types: ", !ismissing(sim.point_types) ? size(sim.point_types) : missing)
    println("  Electric field: ", !ismissing(sim.electric_field) ? size(sim.electric_field) : missing)
    println("  Weighting potentials: ")
    for contact in sim.detector.contacts
        print("    Contact $(contact.id): ")
        println(!ismissing(sim.weighting_potentials[contact.id]) ? size(sim.weighting_potentials[contact.id]) : missing)
    end
    println("  Electron drift field: ", !ismissing(sim.electron_drift_field) ? size(sim.electron_drift_field) : missing)
    println("  Hole drift field: ", !ismissing(sim.hole_drift_field) ? size(sim.hole_drift_field) : missing)
end

function print(io::IO, sim::Simulation{T}) where {T <: SSDFloat}
    print(io, "Simulation{$T} - ", "$(sim.detector.name)")
end

function show(io::IO, sim::Simulation{T}) where {T <: SSDFloat} println(io, sim) end

function show(io::IO, ::MIME"text/plain", sim::Simulation{T}) where {T <: SSDFloat}
    show(io, sim)
end


function Simulation{T}(parsed_dict::Dict)::Simulation{T} where {T <: SSDFloat}
    CS::CoordinateSystemType = Cartesian
    if haskey(parsed_dict, "grid")
        if isa(parsed_dict["grid"], Dict)
            CS = if parsed_dict["grid"]["coordinates"] == "cartesian" 
                Cartesian
            elseif parsed_dict["grid"]["coordinates"]  == "cylindrical"
                Cylindrical
            else
                @assert "`grid` in config file needs `coordinates` that are either `cartesian` or `cylindrical`"
            end
        elseif isa(parsed_dict["grid"], String)
            CS = if parsed_dict["grid"] == "cartesian" 
                Cartesian
            elseif parsed_dict["grid"] == "cylindrical"
                Cylindrical
            else
                @assert "`grid` type in config file needs to be either `cartesian` or `cylindrical`"
            end
        end
    end
    sim::Simulation{T,CS} = Simulation{T,CS}()
    sim.config_dict = parsed_dict
    sim.input_units = construct_units(parsed_dict)
    sim.medium = material_properties[materials[haskey(parsed_dict, "medium") ? parsed_dict["medium"] : "vacuum"]]
    sim.detector = SolidStateDetector{T}(parsed_dict, sim.input_units) 
    sim.world = if haskey(parsed_dict, "grid") && isa(parsed_dict["grid"], Dict)
            World(T, parsed_dict["grid"], sim.input_units)
        else let ssd = sim.detector 
            world_limits = get_world_limits_from_objects(CS, ssd.semiconductor, ssd.contacts, ssd.passives)
            World(CS, world_limits)
        end
    end
    sim.weighting_potentials = Missing[ missing for i in 1:length(sim.detector.contacts)]
    return sim
end

function Simulation{T}(config_file::AbstractString)::Simulation{T} where{T <: SSDFloat}
    parsed_dict = parse_config_file(config_file)
    return Simulation{T}( parsed_dict )
end
function Simulation(config_file::AbstractString)::Simulation{Float32}
    return Simulation{Float32}( config_file )
end

# Functions

function Grid(  sim::Simulation{T, Cylindrical};
                init_grid_size::Union{Missing, NTuple{3, Int}} = missing,
                init_grid_spacing::Union{Missing, Tuple{<:Real,<:Real,<:Real}} = missing,
                for_weighting_potential::Bool = false,
                min_n_ticks::Int = 10,
                full_2π::Bool = false)::CylindricalGrid{T} where {T}
    if ismissing(init_grid_size)
        world_diffs = [(getproperty.(sim.world.intervals, :right) .- getproperty.(sim.world.intervals, :left))...]
        world_diffs[2] = world_diffs[2] * 0.3 * sim.world.intervals[1].right # in radiance
        inds::Vector{Int} = sortperm([world_diffs...])
        ratio::T = min_n_ticks * if world_diffs[inds[1]] > 0
            inv(world_diffs[inds[1]])
        elseif world_diffs[inds[2]] > 0
            inv(world_diffs[inds[2]])
        elseif world_diffs[inds[3]] > 0
            inv(world_diffs[inds[3]])
        else
            error("This should not happen... World has no dimension")
        end
        init_grid_size_1::Int = convert(Int, round(ratio * world_diffs[1], RoundUp))
        init_grid_size_2::Int = convert(Int, round(ratio * world_diffs[2], RoundUp))
        init_grid_size_3::Int = convert(Int, round(ratio * world_diffs[3], RoundUp))
        init_grid_size::NTuple{3, Int} = NTuple{3, T}( [init_grid_size_1, init_grid_size_2, init_grid_size_3] )
    end

    init_grid_spacing, use_spacing::Bool = if !ismissing(init_grid_spacing)
        T.(init_grid_spacing), true
    else
        missing, false
    end
    
    samples::Vector{CylindricalPoint{T}} = sample(sim.detector, Cylindrical)
   
    important_r_points::Vector{T} = map(p -> p.r, samples)
    important_φ_points::Vector{T} = map(p -> p.φ, samples)
    important_z_points::Vector{T} = map(p -> p.z, samples)

    push!(important_r_points, sim.world.intervals[1].left)
    push!(important_r_points, sim.world.intervals[1].right)
    important_r_points = unique!(sort!(geom_round.(important_r_points)))
    push!(important_z_points, sim.world.intervals[3].left)
    push!(important_z_points, sim.world.intervals[3].right)
    important_z_points = unique!(sort!(geom_round.(important_z_points)))
    push!(important_φ_points, sim.world.intervals[2].left)
    push!(important_φ_points, sim.world.intervals[2].right)
    important_φ_points = unique!(sort!(geom_round.(important_φ_points)))

    # r
    L, R, BL, BR = get_boundary_types(sim.world.intervals[1])
    int_r = Interval{L, R, T}(sim.world.intervals[1].left, sim.world.intervals[1].right)
    ax_r::DiscreteAxis{T, BL, BR} = if use_spacing
        DiscreteAxis{BL, BR}(int_r, step = init_grid_spacing[1])
    else
        DiscreteAxis{BL, BR}(int_r, length = init_grid_size[1])
    end
    rticks::Vector{T} = merge_axis_ticks_with_important_ticks(ax_r, important_r_points, atol = minimum(diff(ax_r.ticks))/4)
    ax_r = DiscreteAxis{T, BL, BR}(int_r, rticks)


    # φ
    L, R, BL, BR = get_boundary_types(sim.world.intervals[2])
    int_φ = Interval{L, R, T}(sim.world.intervals[2].left, sim.world.intervals[2].right)
    if full_2π == true || (for_weighting_potential && (sim.world.intervals[2].left != sim.world.intervals[2].right))
        L, R, BL, BR = :closed, :open, :periodic, :periodic
        int_φ = Interval{L, R, T}(0, 2π)
    end
    ax_φ = if int_φ.left == int_φ.right
        DiscreteAxis{T, BL, BR}(int_φ, T[int_φ.left])
    else
        if use_spacing
            DiscreteAxis{BL, BR}(int_φ, step = init_grid_spacing[2])
        else
            DiscreteAxis{BL, BR}(int_φ, length = init_grid_size[2])
        end
    end
    if length(ax_φ) > 1
        φticks::Vector{T} = merge_axis_ticks_with_important_ticks(ax_φ, important_φ_points, atol = minimum(diff(ax_φ.ticks))/4)
        ax_φ = typeof(ax_φ)(int_φ, φticks)
    end
    if isodd(length(ax_φ)) && length(ax_φ) > 1 # must be even
        int_φ = ax_φ.interval
        φticks = ax_φ.ticks
        push!(φticks, geom_round((φticks[end] + φticks[end-1]) * 0.5))
        sort!(φticks)
        ax_φ = typeof(ax_φ)(int_φ, φticks) # must be even
    end
    if length(ax_φ) > 1
        @assert iseven(length(ax_φ)) "CylindricalGrid must have even number of points in φ."
    end

    #z
    L, R, BL, BR = get_boundary_types(sim.world.intervals[3])
    int_z = Interval{L, R, T}(sim.world.intervals[3].left, sim.world.intervals[3].right)
    ax_z::DiscreteAxis{T, BL, BR} = if use_spacing
        DiscreteAxis{BL, BR}(int_z, step = init_grid_spacing[3])
    else
        DiscreteAxis{BL, BR}(int_z, length = init_grid_size[3])
    end
    zticks::Vector{T} = merge_axis_ticks_with_important_ticks(ax_z, important_z_points, atol=minimum(diff(ax_z.ticks))/2)
    ax_z = typeof(ax_z)(int_z, zticks)
    if isodd(length(ax_z)) # must be even
        int_z = ax_z.interval
        zticks = ax_z.ticks
        push!(zticks, geom_round((zticks[end] + zticks[end-1]) * 0.5))
        sort!(zticks)
        ax_z = typeof(ax_z)(int_z, zticks) # must be even
    end
    @assert iseven(length(ax_z)) "CylindricalGrid must have even number of points in z."

    return CylindricalGrid{T}( (ax_r, ax_φ, ax_z) )
end


function Grid(  sim::Simulation{T, Cartesian};
                init_grid_size::Union{Missing, NTuple{3, Int}} = missing,
                init_grid_spacing::Union{Missing, Tuple{<:Real,<:Real,<:Real,}} = missing,
                min_n_ticks::Int = 10,
                for_weighting_potential::Bool = false)::CartesianGrid3D{T} where {T}

    if ismissing(init_grid_size)
        world_diffs = [(getproperty.(sim.world.intervals, :right) .- getproperty.(sim.world.intervals, :left))...]
        inds::Vector{Int} = sortperm([world_diffs...])
        ratio::T = min_n_ticks * if world_diffs[inds[1]] > 0
            inv(world_diffs[inds[1]])
        elseif world_diffs[inds[2]] > 0
            inv(world_diffs[inds[2]])
        elseif world_diffs[inds[3]] > 0
            inv(world_diffs[inds[3]])
        else
            error("This should not happen... World has no dimension")
        end
        init_grid_size_1::Int = convert(Int, round(ratio * world_diffs[1], RoundUp))
        init_grid_size_2::Int = convert(Int, round(ratio * world_diffs[2], RoundUp))
        init_grid_size_3::Int = convert(Int, round(ratio * world_diffs[3], RoundUp))
        init_grid_size::NTuple{3, Int} = NTuple{3, T}( [init_grid_size_1, init_grid_size_2, init_grid_size_3] )
    end

    samples::Vector{CartesianPoint{T}} = sample(sim.detector, Cartesian)
    
    important_x_points::Vector{T} = geom_round.(map(p -> p.x, samples))
    important_y_points::Vector{T} = geom_round.(map(p -> p.y, samples))
    important_z_points::Vector{T} = geom_round.(map(p -> p.z, samples))

    init_grid_spacing, use_spacing::Bool = if !ismissing(init_grid_spacing)
        T.(init_grid_spacing), true
    else
        missing, false
    end

    # x
    L, R, BL, BR = get_boundary_types(sim.world.intervals[1])
    int_x = Interval{L, R, T}(sim.world.intervals[1].left, sim.world.intervals[1].right)
    ax_x::DiscreteAxis{T, BL, BR} = if use_spacing
        DiscreteAxis{BL, BR}(int_x, step = init_grid_spacing[1])
    else
        DiscreteAxis{BL, BR}(int_x, length = init_grid_size[1])
    end
    xticks::Vector{T} = merge_axis_ticks_with_important_ticks(ax_x, important_x_points, atol = minimum(diff(ax_x.ticks)) / 2)
    ax_x = typeof(ax_x)(int_x, xticks)
    if isodd(length(ax_x)) # RedBlack dimension must be of even length
        xticks = ax_x.ticks
        push!(xticks, geom_round((xticks[end] + xticks[end-1]) * 0.5))
        sort!(xticks)
        ax_x = DiscreteAxis{T, BL, BR}(int_x, xticks) # must be even
    end
    @assert iseven(length(ax_x)) "CartesianGrid3D must have even number of points in z."

    # y
    L, R, BL, BR = get_boundary_types(sim.world.intervals[2])
    int_y = Interval{L, R, T}(sim.world.intervals[2].left, sim.world.intervals[2].right)
    ax_y::DiscreteAxis{T, BL, BR} = if use_spacing
        DiscreteAxis{BL, BR}(int_y, step = init_grid_spacing[2])
    else
        DiscreteAxis{BL, BR}(int_y, length = init_grid_size[2])
    end
    yticks::Vector{T} = merge_axis_ticks_with_important_ticks(ax_y, important_y_points, atol = minimum(diff(ax_y.ticks)) / 2)
    ax_y = typeof(ax_y)(int_y, yticks)

    # z
    L, R, BL, BR = get_boundary_types(sim.world.intervals[3])
    int_z = Interval{L, R, T}(sim.world.intervals[3].left, sim.world.intervals[3].right)
    ax_z::DiscreteAxis{T, BL, BR} = if use_spacing
        DiscreteAxis{BL, BR}(int_z, step = init_grid_spacing[3])
    else
        DiscreteAxis{BL, BR}(int_z, length = init_grid_size[3])
    end
    zticks::Vector{T} = merge_axis_ticks_with_important_ticks(ax_z, important_z_points, atol = minimum(diff(ax_z.ticks)) / 2)
    ax_z = typeof(ax_z)(int_z, zticks)

    return CartesianGrid3D{T}( (ax_x, ax_y, ax_z) )
end


"""
    function apply_initial_state!(sim::Simulation{T}, ::Type{ElectricPotential}, grid::Grid{T} = Grid(sim))::Nothing

Applies the initial state of the electric potential calculation.
It overwrites `sim.electric_potential`, `sim.q_eff_imp`, `sim.q_eff_fix`, `sim.ϵ` and `sim.point_types`.
"""
function apply_initial_state!(sim::Simulation{T}, ::Type{ElectricPotential}, grid::Grid{T} = Grid(sim))::Nothing where {T <: SSDFloat}
    fssrb::PotentialSimulationSetupRB{T, 3, 4, get_coordinate_system(sim)} =
        PotentialSimulationSetupRB(sim.detector, grid, sim.medium);

    sim.q_eff_imp = EffectiveChargeDensity(EffectiveChargeDensityArray(fssrb), grid)
    sim.q_eff_fix = EffectiveChargeDensity(FixedEffectiveChargeDensityArray(fssrb), grid)
    sim.ϵ_r = DielectricDistribution(DielektrikumDistributionArray(fssrb), grid)
    sim.point_types = PointTypes(PointTypeArray(fssrb), grid)
    sim.electric_potential = ElectricPotential(ElectricPotentialArray(fssrb), grid)
    nothing
end

"""
    function apply_initial_state!(sim::Simulation{T}, ::Type{WeightingPotential}, contact_id::Int, grid::Grid{T} = Grid(sim))::Nothing

Applies the initial state of the weighting potential calculation for the contact with the id `contact_id`.
It overwrites `sim.weighting_potentials[contact_id]`.
"""
function apply_initial_state!(sim::Simulation{T}, ::Type{WeightingPotential}, contact_id::Int, grid::Grid{T} = Grid(sim))::Nothing where {T <: SSDFloat}
    fssrb::PotentialSimulationSetupRB{T, 3, 4, get_coordinate_system(sim)} =
        PotentialSimulationSetupRB(sim.detector, grid, sim.medium, weighting_potential_contact_id = contact_id);

    sim.weighting_potentials[contact_id] = WeightingPotential(ElectricPotentialArray(fssrb), grid)
    nothing
end


"""
    function update_till_convergence!( sim::Simulation{T} ::Type{ElectricPotential}, convergence_limit::Real; kwargs...)::T

Takes the current state of `sim.electric_potential` and updates it until it has converged.
"""
function update_till_convergence!( sim::Simulation{T,CS},
                                   ::Type{ElectricPotential},
                                   convergence_limit::Real = 1e-7;
                                   n_iterations_between_checks::Int = 500,
                                   max_n_iterations::Int = -1,
                                   depletion_handling::Bool = false,
                                   use_nthreads::Int = Base.Threads.nthreads(),
                                   sor_consts::Union{Missing, T, NTuple{2, T}} = missing
                                    )::T where {T <: SSDFloat, CS <: AbstractCoordinateSystem}
    if ismissing(sor_consts)
        sor_consts = CS == Cylindrical ? (T(1.4), T(1.85)) : T(1.4)
    elseif length(sor_consts) == 1 && CS == Cylindrical
        sor_consts = (T(sor_consts), T(sor_consts))
    elseif length(sor_consts) > 1 && CS == Cartesian
        sor_consts = T(sor_consts[1])
    end
    only_2d::Bool = length(sim.electric_potential.grid.axes[2]) == 1

    fssrb::PotentialSimulationSetupRB{T, 3, 4, get_coordinate_system(sim.electric_potential.grid)} =
        PotentialSimulationSetupRB(sim.detector, sim.electric_potential.grid, sim.medium, sim.electric_potential.data, sor_consts = T.(sor_consts))

    cf::T = _update_till_convergence!( fssrb, T(convergence_limit);
                                       only2d = Val{only_2d}(),
                                       depletion_handling = Val{depletion_handling}(),
                                       is_weighting_potential = Val{false}(),
                                       use_nthreads = use_nthreads,
                                       n_iterations_between_checks = n_iterations_between_checks,
                                       max_n_iterations = max_n_iterations )

    grid::Grid = Grid(fssrb)
    sim.q_eff_imp = EffectiveChargeDensity(EffectiveChargeDensityArray(fssrb), grid)
    sim.q_eff_fix = EffectiveChargeDensity(FixedEffectiveChargeDensityArray(fssrb), grid)
    sim.ϵ_r = DielectricDistribution(DielektrikumDistributionArray(fssrb), grid)
    sim.electric_potential = ElectricPotential(ElectricPotentialArray(fssrb), grid)
    sim.point_types = PointTypes(PointTypeArray(fssrb), grid)

    if depletion_handling == true
        update_again::Bool = false # With SOR-Constant = 1
        @inbounds for i in eachindex(sim.electric_potential.data)
            if sim.electric_potential.data[i] < fssrb.minimum_applied_potential # p-type
                sim.electric_potential.data[i] = fssrb.minimum_applied_potential
                update_again = true
            end
        end
        if update_again == false
            @inbounds for i in eachindex(sim.electric_potential.data)
                if sim.electric_potential.data[i] > fssrb.maximum_applied_potential # n-type
                    sim.electric_potential.data[i] = fssrb.maximum_applied_potential
                    update_again = true
                end
            end
        end
        if update_again
            fssrb.sor_const[:] .= T(1)
            cf = _update_till_convergence!( fssrb, T(convergence_limit);
                                            only2d = Val{only_2d}(),
                                            depletion_handling = Val{depletion_handling}(),
                                            is_weighting_potential = Val{false}(),
                                            use_nthreads = use_nthreads,
                                            n_iterations_between_checks = n_iterations_between_checks,
                                            max_n_iterations = max_n_iterations )
            sim.electric_potential = ElectricPotential(ElectricPotentialArray(fssrb), grid)

            @inbounds for i in eachindex(sim.electric_potential.data)
                if sim.electric_potential.data[i] < fssrb.minimum_applied_potential # p-type
                    @warn """Detector seems not to be fully depleted at a bias voltage of $(fssrb.bias_voltage) V.
                        At least one grid point has a smaller potential value ($(sim.electric_potential.data[i]) V)
                        than the minimum applied potential ($(fssrb.minimum_applied_potential) V). This should not be.
                        However, small overshoots could be due to numerical precision."""
                    break
                end
                if sim.electric_potential.data[i] > fssrb.maximum_applied_potential # n-type
                    @warn """Detector seems not to be not fully depleted at a bias voltage of $(fssrb.bias_voltage) V.
                        At least one grid point has a higher potential value ($(sim.electric_potential.data[i]) V)
                        than the maximum applied potential ($(fssrb.maximum_applied_potential) V). This should not be.
                        However, small overshoots could be due to numerical precision."""
                    break
                end
            end
        end
    end

    cf
end

"""
    function update_till_convergence!( sim::Simulation{T} ::Type{WeightingPotential}, contact_id::Int, convergence_limit::Real; kwargs...)::T

Takes the current state of `sim.weighting_potentials[contact_id]` and updates it until it has converged.
"""
function update_till_convergence!( sim::Simulation{T, CS},
                                   ::Type{WeightingPotential},
                                   contact_id::Int,
                                   convergence_limit::Real;
                                   n_iterations_between_checks::Int = 500,
                                   max_n_iterations::Int = -1,
                                   depletion_handling::Bool = false,
                                   use_nthreads::Int = Base.Threads.nthreads(),
                                   sor_consts::Union{Missing, T, NTuple{2, T}} = missing
                                    )::T where {T <: SSDFloat, CS <: AbstractCoordinateSystem}
    if ismissing(sor_consts)
        sor_consts = CS == Cylindrical ? (T(1.4), T(1.85)) : T(1.4)
    elseif length(sor_consts) == 1 && CS == Cylindrical
        sor_consts = (T(sor_consts), T(sor_consts))
    elseif length(sor_consts) > 1 && CS == Cartesian
        sor_consts = T(sor_consts[1])
    end

    only_2d::Bool = length(sim.weighting_potentials[contact_id].grid.axes[2]) == 1
    fssrb::PotentialSimulationSetupRB{T, 3, 4, get_coordinate_system(sim.weighting_potentials[contact_id].grid)} =
        PotentialSimulationSetupRB(sim.detector, sim.weighting_potentials[contact_id].grid, sim.medium, sim.weighting_potentials[contact_id].data,
                sor_consts = T.(sor_consts), weighting_potential_contact_id = contact_id)

    cf::T = _update_till_convergence!( fssrb, T(convergence_limit);
                                       only2d = Val{only_2d}(),
                                       depletion_handling = Val{depletion_handling}(),
                                       is_weighting_potential = Val{true}(),
                                       use_nthreads = use_nthreads,
                                       n_iterations_between_checks = n_iterations_between_checks,
                                       max_n_iterations = max_n_iterations )

    sim.weighting_potentials[contact_id] = WeightingPotential(ElectricPotentialArray(fssrb), sim.weighting_potentials[contact_id].grid)

    cf
end

"""
    function refine!(sim::Simulation{T}, ::Type{ElectricPotential}, max_diffs::Tuple{<:Real,<:Real,<:Real}, minimum_distances::Tuple{<:Real,<:Real,<:Real})

Takes the current state of `sim.electric_potential` and refines it with respect to the input arguments
`max_diffs` and `minimum_distances`.
"""
function refine!(sim::Simulation{T}, ::Type{ElectricPotential},
                    max_diffs::Tuple{<:Real,<:Real,<:Real} = (T(0), T(0), T(0)),
                    minimum_distances::Tuple{<:Real,<:Real,<:Real} = (T(0), T(0), T(0));
                    update_other_fields::Bool = false) where {T <: SSDFloat}
    sim.electric_potential = refine(sim.electric_potential, max_diffs, minimum_distances)

    if update_other_fields
        fssrb::PotentialSimulationSetupRB{T, 3, 4, get_coordinate_system(sim.electric_potential.grid)} =
            PotentialSimulationSetupRB(sim.detector, sim.electric_potential.grid, sim.medium, sim.electric_potential.data)

        sim.q_eff_imp = EffectiveChargeDensity(EffectiveChargeDensityArray(fssrb), sim.electric_potential.grid)
        sim.q_eff_fix = EffectiveChargeDensity(FixedEffectiveChargeDensityArray(fssrb), sim.electric_potential.grid)
        sim.ϵ_r = DielectricDistribution(DielektrikumDistributionArray(fssrb), sim.electric_potential.grid)
        sim.point_types = PointTypes(PointTypeArray(fssrb), sim.electric_potential.grid)
    end
    nothing
end
"""
    function refine!(sim::Simulation{T}, ::Type{WeightingPotential}, max_diffs::Tuple{<:Real,<:Real,<:Real}, minimum_distances::Tuple{<:Real,<:Real,<:Real})

Takes the current state of `sim.weighting_potentials[contact_id]` and refines it with respect to the input arguments
`max_diffs` and `minimum_distances`.
"""
function refine!(sim::Simulation{T}, ::Type{WeightingPotential}, contact_id::Int,
                    max_diffs::Tuple{<:Real,<:Real,<:Real} = (T(0), T(0), T(0)),
                    minimum_distances::Tuple{<:Real,<:Real,<:Real} = (T(0), T(0), T(0))) where {T <: SSDFloat}
    sim.weighting_potentials[contact_id] = refine(weighting_potentials[contact_id], max_diffs, minimum_distances)
    nothing
end


function _calculate_potential!( sim::Simulation{T, CS}, potential_type::UnionAll, contact_id::Union{Missing, Int} = missing;
        init_grid_size::Union{Missing, NTuple{3, Int}} = missing,
        init_grid_spacing::Union{Missing, Tuple{<:Real,<:Real,<:Real}} = missing,
        grid::Union{Missing, Grid{T}} = missing,
        convergence_limit::Real = 1e-7,
        max_refinements::Int = 3,
        refinement_limits::Union{Missing, Tuple{<:Real,<:Real,<:Real}} = missing,
        min_grid_spacing::Union{Missing, Tuple{<:Real,<:Real,<:Real}} = missing,
        depletion_handling::Bool = false,
        use_nthreads::Int = Base.Threads.nthreads(),
        sor_consts::Union{Missing, <:Real, Tuple{<:Real,<:Real}} = missing,
        max_n_iterations::Int = 50000,
        verbose::Bool = true,
    )::Nothing where {T <: SSDFloat, CS <: AbstractCoordinateSystem}

    begin # preperations
        convergence_limit::T = T(convergence_limit)
        isEP::Bool = potential_type == ElectricPotential
        isWP::Bool = !isEP
        if isWP depletion_handling = false end
        if ismissing(grid)
            grid = Grid(sim, init_grid_size = init_grid_size, init_grid_spacing = init_grid_spacing, for_weighting_potential = isWP)
        end
        if ismissing(sor_consts)
            sor_consts = CS == Cylindrical ? (T(1.4), T(1.85)) : T(1.4)
        elseif length(sor_consts) == 1 && CS == Cylindrical
            sor_consts = (T(sor_consts), T(sor_consts))
        elseif length(sor_consts) > 1 && CS == Cartesian
            sor_consts = T(sor_consts[1])
        end
        if ismissing(refinement_limits)
            refinement_limits::NTuple{3, T} = CS == Cylindrical ? (T(1e-5), T(1e-5) / (0.25 * grid.axes[1][end]), T(1e-5)) : (T(1e-5), T(1e-5), T(1e-5))
        end
        if ismissing(min_grid_spacing)
            min_grid_spacing::NTuple{3, T} = CS == Cylindrical ? (T(1e-5), T(1e-5) / (0.25 * grid.axes[1][end]), T(1e-5)) : (T(1e-5), T(1e-5), T(1e-5))
        end
        n_iterations_between_checks::Int = div(10^7, length(grid))
        if n_iterations_between_checks < 20 n_iterations_between_checks = 20 end
        if use_nthreads > Base.Threads.nthreads()
            use_nthreads = Base.Threads.nthreads();
            @warn "`use_nthreads` was set to `1`. The environment variable `JULIA_NUM_THREADS` must be set appropriately before the julia session is started."
        end
        refine::Bool = max_refinements > 0 ? true : false
        only_2d::Bool = length(grid.axes[2]) == 1 ? true : false
        if CS == Cylindrical
            cyclic::T = grid.axes[2].interval.right - grid.axes[2].interval.left
            n_φ_sym::Int = only_2d ? 1 : round(Int, round(T(2π) / cyclic, digits = 3))
            n_φ_sym_info_txt = if only_2d
                "φ symmetry: Detector is φ-symmetric -> 2D computation."
            else
                "φ symmetry: calculating just 1/$(n_φ_sym) in φ of the detector."
            end
        end
        contact_potentials::Vector{T} = [contact.potential for contact in sim.detector.contacts]
        bias_voltage::T = (length(contact_potentials) > 0) ? (maximum(contact_potentials) - minimum(contact_potentials)) : T(0)
        if isWP bias_voltage = T(1) end
        if verbose
            println(
                "$(isEP ? "Electric" : "Weighting") Potential Calculation\n",
                if isEP "Bias voltage: $(bias_voltage) V\n" else "" end,
                if CS == Cylindrical "$n_φ_sym_info_txt\n" else "" end,
                "Precision: $T\n",
                "Convergence limit: $convergence_limit => $(round(abs(bias_voltage * convergence_limit), sigdigits=2)) V\n",
                "Threads: $use_nthreads\n",
                "Coordinate system: $(CS)\n",
                "Initial grid dimension: $(size(grid))\n",
                "Refine? -> $refine\n",
                "Refinement parameters:\n",
                "\tmaximum number of refinements: $(max_refinements)\n",
                "\tminimum grid spacing:\n",
                "\t\tr: $(min_grid_spacing[1]) m\n",
                "\t\tφ: $(min_grid_spacing[2]) rad\n",
                "\t\tz: $(min_grid_spacing[3]) m\n",
                "\tRefinement limits:\n",
                "\t\tr: $(refinement_limits[1]) -> $(round(abs(bias_voltage * refinement_limits[1]), sigdigits=2)) V\n",
                "\t\tφ: $(refinement_limits[2]) -> $(round(abs(bias_voltage * refinement_limits[2]), sigdigits=2)) V\n",
                "\t\tz: $(refinement_limits[3]) -> $(round(abs(bias_voltage * refinement_limits[3]), sigdigits=2)) V\n",
                ""
            )
        end
    end
    if isEP
        apply_initial_state!(sim, potential_type, grid)
        update_till_convergence!( sim, potential_type, convergence_limit,
                                  n_iterations_between_checks = n_iterations_between_checks,
                                  max_n_iterations = max_n_iterations,
                                  depletion_handling = depletion_handling,
                                  use_nthreads = use_nthreads,
                                  sor_consts = sor_consts )
    else
        apply_initial_state!(sim, potential_type, contact_id, grid)
        update_till_convergence!( sim, potential_type, contact_id, convergence_limit,
                                    n_iterations_between_checks = n_iterations_between_checks,
                                    max_n_iterations = max_n_iterations,
                                    depletion_handling = depletion_handling,
                                    use_nthreads = use_nthreads,
                                    sor_consts = sor_consts )
    end

    refinement_counter::Int = 1
    if refine
        max_diffs::NTuple{3, T} = abs.(refinement_limits .* bias_voltage)
        refine_at_inds::NTuple{3, Vector{Int}} = if isEP
            _get_refinement_inds(sim.electric_potential.data, sim.electric_potential.grid, max_diffs, min_grid_spacing)
        else
            _get_refinement_inds(sim.weighting_potentials[contact_id].data, sim.weighting_potentials[contact_id].grid, max_diffs, min_grid_spacing)
        end
        while any(!isempty, refine_at_inds[1]) && refinement_counter <= max_refinements
            if isEP
                sim.electric_potential = ElectricPotential(add_points_and_interpolate(
                        sim.electric_potential.data, sim.electric_potential.grid, refine_at_inds...)...)
                if verbose @info "New Grid Size = $(size(sim.electric_potential.grid))" end
            else
                sim.weighting_potentials[contact_id] = WeightingPotential(add_points_and_interpolate(
                    sim.weighting_potentials[contact_id].data, sim.weighting_potentials[contact_id].grid, refine_at_inds...)...)
                if verbose @info "New Grid Size = $(size(sim.weighting_potentials[contact_id].grid))" end
            end
            n_iterations_between_checks = div(10^7, length(isEP ? sim.electric_potential.grid : sim.weighting_potentials[contact_id].grid))
            if n_iterations_between_checks < 20 n_iterations_between_checks = 20 end

            if isEP
                update_till_convergence!( sim, potential_type, convergence_limit,
                    n_iterations_between_checks = n_iterations_between_checks,
                    max_n_iterations = max_n_iterations,
                    depletion_handling = depletion_handling,
                    use_nthreads = use_nthreads,
                    sor_consts = sor_consts )
            else
                update_till_convergence!( sim, potential_type, contact_id, convergence_limit,
                                    n_iterations_between_checks = n_iterations_between_checks,
                                    max_n_iterations = max_n_iterations,
                                    depletion_handling = depletion_handling,
                                    use_nthreads = use_nthreads,
                                    sor_consts = sor_consts )
            end

            refine_at_inds = if isEP
                _get_refinement_inds(sim.electric_potential.data, sim.electric_potential.grid, max_diffs, min_grid_spacing)
            else
                _get_refinement_inds(sim.weighting_potentials[contact_id].data, sim.weighting_potentials[contact_id].grid, max_diffs, min_grid_spacing)
            end
            refinement_counter += 1
        end
    end

    nothing
end

"""
    calculate_weighting_potential!(sim::Simulation{T}, contact_id::Int; kwargs...)::Nothing

Compute the weighting potential for the contact with id `contact_id`
for the given Simulation `sim` on an adaptive grid through successive over relaxation.

There are serveral `<keyword arguments>` which can be used to tune the computation:

# Keywords
- `convergence_limit::Real`: `convergence_limit` times the bias voltage sets the convergence limit of the relaxation.
    The convergence value is the absolute maximum difference of the potential between two iterations of all grid points.
    Default of `convergence_limit` is `2e-6` (times bias voltage).
- `max_refinements::Int`: Number of maximum refinements. Default is `2`. Set it to `0` to switch off refinement.
- `refinement_limits::Tuple{<:Real, <:Real, <:Real}`: Tuple of refinement limits for each dimension
    (in case of cylindrical coordinates the order is `r`, `φ`, `z`).
    A refinement limit (e.g. `refinement_limits[1]`) times the bias voltage of the detector `det` is the
    maximum allowed voltage difference between two neighbouring grid points in the respective dimension.
    When the difference is larger, new points are created inbetween. Default is `[1e-5, 1e-5, 1e-5]`.
- `init_grid_spacing::Tuple{<:Real, <:Real, <:Real}`: Tuple of the initial distances between two grid points for each dimension.
    For normal coordinates the unit is meter. For angular coordinates, the unit is radiance.
    It prevents the refinement to make the grid to fine. Default is `[0.005, 10.0, 0.005]``.
- `min_grid_spacing::Tuple{<:Real, <:Real, <:Real}`: Tuple of the mimimum allowed distance between two grid points for each dimension.
    For normal coordinates the unit is meter. For angular coordinates, the unit is radiance.
    It prevents the refinement to make the grid to fine. Default is [`1e-6`, `1e-6`, `1e-6`].
- `grid::Grid{T, N, S}`: Initial grid used to start the simulation. Default is `Grid(sim, init_grid_spacing=init_grid_spacing)`.
- `depletion_handling::Bool`: Enables the handling of undepleted regions. Default is false.
- `use_nthreads::Int`: Number of threads to use in the computation. Default is `Base.Threads.nthreads()`.
    The environment variable `JULIA_NUM_THREADS` must be set appropriately before the Julia session was
    started (e.g. `export JULIA_NUM_THREADS=8` in case of bash).
- `sor_consts::Union{<:Real, NTuple{2, <:Real}}`: Two element tuple in case of cylindrical coordinates.
    First element contains the SOR constant for `r` = 0.
    Second contains the constant at the outer most grid point in `r`. A linear scaling is applied in between.
    First element should be smaller than the second one and both should be ∈ [1.0, 2.0]. Default is [1.4, 1.85].
    In case of cartesian coordinates only one value is taken.
- `max_n_iterations::Int`: Set the maximum number of iterations which are performed after each grid refinement.
    Default is `10000`. If set to `-1` there will be no limit.
- `verbose::Bool=true`: Boolean whether info output is produced or not.
"""
function calculate_weighting_potential!(sim::Simulation{T}, contact_id::Int, args...; n_points_in_φ::Union{Missing, Int} = missing, kwargs...)::Nothing where {T <: SSDFloat}
    # S = get_coordinate_system(sim)
    # periodicity::T = get_periodicity(sim.world.intervals[2])
    # if S == Cylindrical && periodicity == T(0)
    #     if ismissing(n_points_in_φ)
    #         @info "\tIn weighing potential calculation: Keyword `n_points_in_φ` not set.\n\t\tDefault is `n_points_in_φ = 36`. 2D field will be extended to 36 points in φ."
    #         n_points_in_φ = 36
    #     else
    #         if !(n_points_in_φ > 1 && iseven(n_points_in_φ))
    #             @info "\tIn weighing potential calculation: Keyword `n_points_in_φ` is $(n_points_in_φ) but must be even and larger than 1.\n\t\t`n_points_in_φ` is now set to 36. 2D field will be extended to 36 points in φ."
    #             n_points_in_φ = 36
    #         end
    #     end
    # end
    # wps = calculate_weighting_potential(sim.detector, contact_id, args...; kwargs...)
    # if S == Cylindrical && size(wps.potential, 2) == 1 && !ismissing(n_points_in_φ)
    #     wp = WeightingPotential(wps, n_points_in_φ = n_points_in_φ)
    # else
    #     wp = WeightingPotential(wps)
    # end
    _calculate_potential!(sim, WeightingPotential, contact_id, args...; kwargs...)
    nothing
end

"""
    calculate_electric_potential!(sim::Simulation{T}; kwargs...)::Nothing


Compute the electric potential for the given Simulation `sim` on an adaptive grid
through successive over relaxation.

There are serveral `<keyword arguments>` which can be used to tune the computation:

# Keywords
- `convergence_limit::Real`: `convergence_limit` times the bias voltage sets the convergence limit of the relaxation.
    The convergence value is the absolute maximum difference of the potential between two iterations of all grid points.
    Default of `convergence_limit` is `2e-6` (times bias voltage).
- `max_refinements::Int`: Number of maximum refinements. Default is `2`. Set it to `0` to switch off refinement.
- `refinement_limits::Tuple{<:Real, <:Real, <:Real}`: Tuple of refinement limits for each dimension
    (in case of cylindrical coordinates the order is `r`, `φ`, `z`).
    A refinement limit (e.g. `refinement_limits[1]`) times the bias voltage of the detector `det` is the
    maximum allowed voltage difference between two neighbouring grid points in the respective dimension.
    When the difference is larger, new points are created inbetween. Default is `[1e-5, 1e-5, 1e-5]`.
- `init_grid_spacing::Tuple{<:Real, <:Real, <:Real}`: Tuple of the initial distances between two grid points for each dimension.
    For normal coordinates the unit is meter. For angular coordinates, the unit is radiance.
    It prevents the refinement to make the grid to fine. Default is `[0.005, 10.0, 0.005]``.
- `min_grid_spacing::Tuple{<:Real, <:Real, <:Real}`: Tuple of the mimimum allowed distance between two grid points for each dimension.
    For normal coordinates the unit is meter. For angular coordinates, the unit is radiance.
    It prevents the refinement to make the grid to fine. Default is [`1e-6`, `1e-6`, `1e-6`].
- `grid::Grid{T, N, S}`: Initial grid used to start the simulation. Default is `Grid(sim, init_grid_spacing=init_grid_spacing)`.
- `depletion_handling::Bool`: Enables the handling of undepleted regions. Default is false.
- `use_nthreads::Int`: Number of threads to use in the computation. Default is `Base.Threads.nthreads()`.
    The environment variable `JULIA_NUM_THREADS` must be set appropriately before the Julia session was
    started (e.g. `export JULIA_NUM_THREADS=8` in case of bash).
- `sor_consts::Union{<:Real, NTuple{2, <:Real}}`: Two element tuple in case of cylindrical coordinates.
    First element contains the SOR constant for `r` = 0.
    Second contains the constant at the outer most grid point in `r`. A linear scaling is applied in between.
    First element should be smaller than the second one and both should be ∈ [1.0, 2.0]. Default is [1.4, 1.85].
    In case of cartesian coordinates only one value is taken.
- `max_n_iterations::Int`: Set the maximum number of iterations which are performed after each grid refinement.
    Default is `10000`. If set to `-1` there will be no limit.
- `verbose::Bool=true`: Boolean whether info output is produced or not.
"""
function calculate_electric_potential!(sim::Simulation{T}, args...; kwargs...)::Nothing where {T <: SSDFloat}
    _calculate_potential!(sim, ElectricPotential, args...; kwargs...)
    nothing
end

"""
    calculate_electric_field!(sim::Simulation{T}, args...; n_points_in_φ::Union{Missing, Int} = missing, kwargs...)::Nothing

ToDo...
"""
function calculate_electric_field!(sim::Simulation{T, CS}, args...; n_points_in_φ::Union{Missing, Int} = missing, kwargs...)::Nothing where {T <: SSDFloat, CS}
    periodicity::T = get_periodicity(sim.world.intervals[2])
    e_pot, point_types = if CS == Cylindrical && periodicity == T(0) # 2D, only one point in φ
        if ismissing(n_points_in_φ)
            @info "\tIn electric field calculation: Keyword `n_points_in_φ` not set.\n\t\tDefault is `n_points_in_φ = 36`. 2D field will be extended to 36 points in φ."
            n_points_in_φ = 36
        else
            if !(n_points_in_φ > 1 && iseven(n_points_in_φ))
                @info "\tIn electric field calculation: Keyword `n_points_in_φ` is $(n_points_in_φ) but must be even and larger than 1.\n\t\t`n_points_in_φ` is now set to 36. 2D field will be extended to 36 points in φ."
                n_points_in_φ = 36
            end
        end
        get_2π_potential(sim.electric_potential, n_points_in_φ = n_points_in_φ),
        get_2π_potential(sim.point_types,  n_points_in_φ = n_points_in_φ);
    elseif CS == Cylindrical
        get_2π_potential(sim.electric_potential),
        get_2π_potential(sim.point_types)
    else
        sim.electric_potential,
        sim.point_types
    end
    sim.electric_field = get_electric_field_from_potential(e_pot, point_types);
    nothing
end

function set_charge_drift_model!(sim::Simulation{T}, charge_drift_model::AbstractChargeDriftModel{T})::Nothing where {T <: SSDFloat}
    sim.detector.semiconductor.charge_drift_model = charge_drift_model
    nothing
end

function calculate_drift_fields!(sim::Simulation{T};
    use_nthreads::Int = Base.Threads.nthreads())::Nothing where {T <: SSDFloat}
    sim.electron_drift_field = ElectricField(get_electron_drift_field(sim.electric_field.data, sim.detector.semiconductor.charge_drift_model, use_nthreads = use_nthreads), sim.electric_field.grid)
    sim.hole_drift_field = ElectricField(get_hole_drift_field(sim.electric_field.data, sim.detector.semiconductor.charge_drift_model, use_nthreads = use_nthreads), sim.electric_field.grid)
    nothing
end
@deprecate apply_charge_drift_model!(args...; kwargs...) calculate_drift_fields!(args...; kwargs...)

function get_interpolated_drift_field(ef::ElectricField)
    get_interpolated_drift_field(ef.data, ef.grid)
end

function drift_charges( sim::Simulation{T}, starting_positions::Vector{CartesianPoint{T}};
                        Δt::RealQuantity = 5u"ns", max_nsteps::Int = 1000, verbose::Bool = true )::Vector{EHDriftPath{T}} where {T <: SSDFloat}
    return _drift_charges(   sim.detector, sim.point_types.grid, sim.point_types, starting_positions,
                             get_interpolated_drift_field(sim.electron_drift_field), get_interpolated_drift_field(sim.hole_drift_field),
                             Δt, max_nsteps = max_nsteps, verbose = verbose)::Vector{EHDriftPath{T}}
end

function get_signal(sim::Simulation{T, CS}, drift_paths::Vector{EHDriftPath{T}}, energy_depositions::Vector{T}, contact_id::Int; Δt::TT = T(5) * u"ns") where {T <: SSDFloat, CS, TT}
    dt::T = to_internal_units(internal_time_unit, Δt)
    wp::Interpolations.Extrapolation{T, 3} = interpolated_scalarfield(sim.weighting_potentials[contact_id])
    timestamps = _common_timestamps( drift_paths, dt )
    signal::Vector{T} = zeros(T, length(timestamps))
    add_signal!(signal, timestamps, drift_paths, energy_depositions, wp, CS)
    return RDWaveform( range(zero(T) * unit(Δt), step = T(ustrip(Δt)) * unit(Δt), length = length(signal)), signal )
end

"""
    function simulate!(sim::Simulation{T};  max_refinements::Int = 1, verbose::Bool = false,
                                        depletion_handling::Bool = false, convergence_limit::Real = 1e-5 ) where {T <: SSDFloat}

ToDo...
"""
function simulate!(sim::Simulation{T};  max_refinements::Int = 1, verbose::Bool = false,
                                        depletion_handling::Bool = false, convergence_limit::Real = 1e-7 ) where {T <: SSDFloat}
    calculate_electric_potential!(  sim,
                                    max_refinements = max_refinements,
                                    verbose = verbose,
                                    init_grid_size = (10,10,10),
                                    depletion_handling = depletion_handling,
                                    convergence_limit = convergence_limit )
    for contact in sim.detector.contacts
        calculate_weighting_potential!(sim, contact.id, max_refinements = max_refinements,
                init_grid_size = (10,10,10),
                verbose = verbose, convergence_limit = convergence_limit)
    end
    calculate_electric_field!(sim)
    calculate_drift_fields!(sim)
    @info "Detector simulation done"
end

function _get_abs_bias_voltage(det::SolidStateDetector{T}) where {T <: SSDFloat}
    potentials::Vector{T} = map(c -> c.potential, det.contacts)
    return (maximum(potentials) - minimum(potentials)) * u"V"
end




calculate_stored_energy(sim::Simulation) =
    calculate_stored_energy(sim.electric_field, sim.ϵ_r)

function calculate_stored_energy(ef::ElectricField{T,3,S}, ϵ::DielectricDistribution{T,3,S}) where {T <: SSDFloat, S}
    W::T = 0

    cylindric::Bool = S == Cylindrical
    cartesian::Bool = !cylindric
    r0_handling::Bool = typeof(ef.grid.axes[1]).parameters[2] == :r0

    ax1::Vector{T} = collect(ef.grid[1])
    ax2::Vector{T} = collect(ef.grid[2])
    ax3::Vector{T} = collect(ef.grid[3])
    mp1::Vector{T} = midpoints(get_extended_ticks(ef.grid[1]))
    mp2::Vector{T} = midpoints(get_extended_ticks(ef.grid[2]))
    mp3::Vector{T} = midpoints(get_extended_ticks(ef.grid[3]))
    Δax1::Vector{T} = diff(ax1)
    Δax2::Vector{T} = diff(ax2)
    Δax3::Vector{T} = diff(ax3)
    Δmp1::Vector{T} = diff(mp1)
    Δmp2::Vector{T} = diff(mp2)
    Δmp3::Vector{T} = diff(mp3)

    w1r::Vector{T} = inv.(Δmp1) .* (mp1[2:end] .- ax1)
    w1l::Vector{T} = inv.(Δmp1) .* (ax1 - mp1[1:end-1])
    w2r::Vector{T} = inv.(Δmp2) .* (mp2[2:end] .- ax2)
    w2l::Vector{T} = inv.(Δmp2) .* (ax2 - mp2[1:end-1])
    w3r::Vector{T} = inv.(Δmp3) .* (mp3[2:end] .- ax3)
    w3l::Vector{T} = inv.(Δmp3) .* (ax3 - mp3[1:end-1])

    if cylindric
        mp1[1] = 0
        mp1[end] = ax1[end]
        Δmp1 = ((mp1[2:end].^2) .- (mp1[1:end-1].^2)) ./ 2
    end
    V::T = 0
    for i3 in 1:size(ϵ, 3)-1
        _Δmp3::T = Δmp3[i3]
        if (i3 == 1 || i3 == size(ϵ, 3)-1) _Δmp3 /= 2 end
        for i2 in 1:size(ϵ, 2)-1
            _Δmp2::T = Δmp2[i2]
            if (cartesian && (i2 == 1 || i2 == size(ϵ, 2)-1)) _Δmp2 /= 2 end
            for i1 in 1:size(ϵ, 1)-1
                _Δmp1::T = Δmp1[i1]
                if (cartesian && (i1 == 1 || i1 == size(ϵ, 1)-1)) _Δmp1 /= 2 end
                ev::SArray{Tuple{3},Float32,1,3} = ef.data[i1, i2, i3]
                dV::T = _Δmp3 * _Δmp2 * _Δmp1
                _ϵ::T = sum([
                    ϵ[i1, i2, i3]             * w1l[i1] * w2l[i2] * w3l[i3],
                    ϵ[i1 + 1, i2, i3]         * w1r[i1] * w2l[i2] * w3l[i3],
                    ϵ[i1, i2 + 1, i3]         * w1l[i1] * w2r[i2] * w3l[i3],
                    ϵ[i1, i2, i3 + 1]         * w1l[i1] * w2l[i2] * w3r[i3],
                    ϵ[i1 + 1, i2 + 1, i3]     * w1r[i1] * w2r[i2] * w3l[i3],
                    ϵ[i1 + 1, i2, i3 + 1]     * w1r[i1] * w2l[i2] * w3r[i3],
                    ϵ[i1, i2 + 1, i3 + 1]     * w1l[i1] * w2r[i2] * w3r[i3],
                    ϵ[i1 + 1, i2 + 1, i3 + 1] * w1r[i1] * w2r[i2] * w3r[i3]
                ])
                V += dV
                W += sum(ev.^2) * dV * _ϵ
            end
        end
    end
    E =  W * ϵ0 / 2 * u"J"
    if cylindric && (size(ϵ, 2) - 1 != size(ef, 2))
        E *= size(ef, 2) / (size(ϵ, 2) - 1)
    end
    return E
end

export calculate_capacitance
"""
    calculate_capacitance(sim::Simulation{T})::T where {T <: SSDFloat}

Calculates the capacitance of an detector in Farad.
"""
function calculate_capacitance(sim::Simulation{T}) where {T <: SSDFloat}
    W = calculate_stored_energy(sim)
    return uconvert(u"pF", 2 * W / (_get_abs_bias_voltage(sim.detector)^2))
end
