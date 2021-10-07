## Electric field strength plotting

@recipe function f( ef::ElectricField{T, 3, Cylindrical};
                    r = missing,
                    φ = missing,
                    z = missing,
                    contours_equal_potential=false,
                    full_det = false ) where {T}

    grid::CylindricalGrid{T} = ef.grid
    cross_section::Symbol, idx::Int, value::T, units::Unitful.Units = get_crosssection_idx_and_value(grid, r, φ, z)
    
    title --> "Electric Field (Magn.) @ $(cross_section) = $(round(units, Float64(uconvert(units,value*(cross_section == :φ ? u"°" : internal_length_unit))), sigdigits=2))"
    colorbar_title --> "Electric Field Strength"
    seriescolor --> :inferno

    ef_magn = norm.(ef)
    ElectricPotential(ef_magn, grid),  cross_section, idx, value, internal_efield_unit, contours_equal_potential, full_det
end

@recipe function f(ef::ElectricField{T, 3, Cartesian};
                    x = missing,
                    y = missing,
                    z = missing,
                    contours_equal_potential = false) where {T <: SSDFloat}

    grid::CartesianGrid3D{T} = ef.grid
    cross_section::Symbol, idx::Int, value::T, units::Unitful.Units = get_crosssection_idx_and_value(grid, x, y, z)
    
    title --> "Electric Field (Magn.) @ $(cross_section) = $(round(units, Float64(uconvert(units,value*internal_length_unit)), sigdigits=2))"
    colorbar_title --> "Electric Field Strength"
    seriescolor --> :inferno

    ef_magn = norm.(ef)
    ElectricPotential(ef_magn, grid), cross_section, idx, value, internal_efield_unit, contours_equal_potential
end


## Electric field line plotting

function get_sample_lines(dim_symbol::Symbol, v::T, grid::CartesianGrid3D{T}, sampling::T)::Vector{ConstructiveSolidGeometry.Line} where {T}
    
    xrange = range(round.(endpoints(grid.x.interval)./sampling, (RoundDown, RoundUp)).*sampling..., step = sampling)
    yrange = range(round.(endpoints(grid.y.interval)./sampling, (RoundDown, RoundUp)).*sampling..., step = sampling)
    zrange = range(round.(endpoints(grid.z.interval)./sampling, (RoundDown, RoundUp)).*sampling..., step = sampling)
    
    return if dim_symbol == :x
        vcat(
            [ConstructiveSolidGeometry.Line(CartesianPoint{T}(v,y,0), CartesianVector{T}(0,0,1)) for y in yrange],
            [ConstructiveSolidGeometry.Line(CartesianPoint{T}(v,0,z), CartesianVector{T}(0,1,0)) for z in zrange]
        )
    elseif dim_symbol == :y
        vcat(
            [ConstructiveSolidGeometry.Line(CartesianPoint{T}(x,v,0), CartesianVector{T}(0,0,1)) for x in xrange],
            [ConstructiveSolidGeometry.Line(CartesianPoint{T}(0,v,z), CartesianVector{T}(1,0,0)) for z in zrange]
        )
        
    else # dim_symbol == :z
        vcat(
            [ConstructiveSolidGeometry.Line(CartesianPoint{T}(x,0,v), CartesianVector{T}(0,1,0)) for x in xrange],
            [ConstructiveSolidGeometry.Line(CartesianPoint{T}(0,y,v), CartesianVector{T}(1,0,0)) for y in yrange]
        )
    end
end


function get_sample_lines(dim_symbol::Symbol, v::T, grid::CylindricalGrid{T}, sampling::T)::Vector{ConstructiveSolidGeometry.Line} where {T}
    
    φsampling = 2π * sampling / (dim_symbol == :r ? v : grid.r.interval.right) # or use sampling together with the unit ?
    rrange = range(round.((-grid.r.interval.right,grid.r.interval.right)./sampling, (RoundDown, RoundUp)).*sampling..., step = sampling)
    φrange = range(round.(endpoints(grid.φ.interval)./φsampling, (RoundDown, RoundUp)).*φsampling..., step = φsampling)
    zrange = range(round.(endpoints(grid.z.interval)./sampling, (RoundDown, RoundUp)).*sampling..., step = sampling)
    
    return if dim_symbol == :φ
        vcat(
            [ConstructiveSolidGeometry.Line(CartesianPoint{T}(r*cos(v),r*sin(v),0), CartesianVector{T}(0,0,1)) for r in rrange],
            [ConstructiveSolidGeometry.Line(CartesianPoint{T}(0,0,z), CartesianVector{T}(cos(v),sin(v),0)) for z in zrange]
        )
    elseif dim_symbol == :z
        [ConstructiveSolidGeometry.Line(CartesianPoint{T}(0,0,v), CartesianVector{T}(cos(φ), sin(φ),0)) for φ in φrange]
    else # dim_symbol == :r
        [ConstructiveSolidGeometry.Line(CartesianPoint{T}(v*cos(φ),v*sin(φ),0), CartesianVector{T}(0,0,1)) for φ in φrange]
    end
end


@userplot Plot_electric_fieldlines
@recipe function f(gdd::Plot_electric_fieldlines; φ = missing, r = missing, x = missing, y = missing, z = missing,
                    max_nsteps = 5000,
                    sampling = 2u"mm", # Specifies in what density the Contacts are sampled to generate equally spaced surface charges. Also see spacing.
                    offset = 0.5u"mm", # should be at least as big as sampling. In doubt sampling can be reduced and the spacing keyword can be used to thin out the lines.
                    spacing = 1, # If, due to fine sampling, too many lines would clutter the plot, the spacing keyword allows to skip some fieldlines. Spacing = 2 means plot every second line. Spacing = 3 every third.
                    full_det = false,
                    skip_contact = 1) # Usually the "core" contact is skipped, and the other contacts are equally sampled for charges to drift
    sim = gdd.args[1]
    S = get_coordinate_system(sim.electric_field.grid)
    T = get_precision_type(sim.detector)
    dim_array::Vector{Union{Missing, RealQuantity}} = (S == Cartesian ? [x, y, z] : [r, φ, z])
    dim_symbols_array = (S == Cartesian ? [:x, :y, :z] : [:r, :φ, :z])

    if isempty(skipmissing(dim_array))
        if S == Cartesian
            dim_number = 1
            dim_symbol = :x
            v::T = sim.electric_field.grid[dim_number][ div(length(sim.electric_field.grid[dim_number]), 2) ]
            x = v
            units::Unitful.Units = u"m"
        else # S == Cylindrical
            v = 0
            φ = 0
            dim_number = 2
            dim_symbol = :φ
            units = u"°"
        end
        dim_array = (S == Cartesian ? [x, y, z] : [r, φ, z])
    elseif sum(ismissing.(dim_array)) == 2
        dim_number = findfirst(x -> !ismissing(x), dim_array)
        dim_symbol = dim_symbols_array[dim_number]
        vtmp = dim_array[dim_number]
        units = vtmp isa Real ? (dim_symbol == :φ ? internal_angle_unit : internal_length_unit) : unit(vtmp)
        v = T(to_internal_units(dim_symbol == :φ && φ isa Real ? deg2rad(vtmp) : vtmp))
    else
        throw(ArgumentError("Only one keyword for a certain dimension is allowed. Please choose one of "*
            (S == Cartesian ? "'x', 'y', 'z'." : "'r', 'φ', 'z'.")))
    end

    show_full_det = full_det && dim_symbol == :φ # The full_det keyword only makes sense for crossections in the xz plane in cylindrical grids

    (dim_symbol != :r && !(dim_symbol == :z && S == Cylindrical) ) ? aspect_ratio --> 1 : nothing
    title --> "Electric Field Lines @ $(dim_symbol) = $(round(units, Float64(uconvert(units, v*(dim_symbol == :φ ? u"rad" : internal_length_unit))), sigdigits=3))" * (show_full_det ? "\n(= $(round(units, Float64(uconvert(units, (v+π)%2π * u"rad")), sigdigits=3)) on left side)" : "")
    xguide --> (S == Cylindrical ? (dim_symbol == :r ? "φ" : "r") : (dim_symbol == :x ? "y" : "x"))
    yguide --> "z"
    unitformat --> :slash
    (S == Cylindrical && dim_symbol == :z) ? xguide :=  "" : nothing
    (S == Cylindrical && dim_symbol == :z) ? yguide :=  "" : nothing

    contacts_to_spawn_charges_for = Contact{T}[c for c in sim.detector.contacts if c.id != skip_contact]
    spawn_positions = CartesianPoint{T}[]
    grid = sim.electric_field.grid
    
    sample_lines = get_sample_lines(dim_symbol, v, grid, T(ustrip(to_internal_units(sampling))))

    for c in contacts_to_spawn_charges_for[:]
        surfs = ConstructiveSolidGeometry.surfaces(c.geometry)
        for l in sample_lines
            for surf in surfs
                pts = ConstructiveSolidGeometry.intersection(surf,l)
                for pt in pts
                    if pt in c
                        pt_in = pt + T(ustrip(to_internal_units(offset))) * normalize(ConstructiveSolidGeometry.normal(surf, pt))
                        if pt_in in sim.detector && !(pt_in in sim.detector.contacts)
                            push!(spawn_positions, pt_in)
                        end
                    end
                end
            end
        end
    end

    el_field_itp     = get_interpolated_drift_field(sim.electric_field.data      , sim.electric_field.grid) #Interpolate the Electric fields, in which the charges will drift. Is passed to the drift function.
    el_field_itp_inv = get_interpolated_drift_field(sim.electric_field.data .* -1, sim.electric_field.grid)

    for (ipos, pos) in enumerate(spawn_positions) # Charge drift and plotting loop. Not optimized for speed, but it doesnt have to be. Uses low level drift function for more contorl.
        for el_field in (el_field_itp, el_field_itp_inv)
            if ((spacing-1)+ipos)%spacing == 0
                path = CartesianPoint{T}[CartesianPoint{T}(0.0,0.0,0.0) for i in 1:max_nsteps]
                _drift_charge!(path, Vector{T}(undef, max_nsteps), sim.detector, sim.point_types, sim.electric_potential.grid, CartesianPoint(pos), T(1e-9), el_field, verbose = false )
                filter!(x->x != CartesianPoint{T}(0.0,0.0,0.0), path)
                @series begin
                    seriescolor --> :white
                    label --> ""
                    x, y = if dim_symbol == :φ
                        map(x -> x[1] * cos(v) + x[2] * sin(v), path), #project to the φ-plane
                        map(x -> x[3], path)
                    elseif dim_symbol == :x
                        map(x -> x[2], path),
                        map(x -> x[3], path)
                    elseif dim_symbol == :y
                        map(x -> x[1], path),
                        map(x -> x[3], path)
                    elseif dim_symbol == :z
                        if S == Cylindrical
                            projection --> :polar
                            path = CylindricalPoint.(path)
                            map(x -> x[2], path),
                            map(x -> x[1], path)
                        else
                            map(x -> x[1], path),
                            map(x -> x[2], path)
                        end
                    elseif dim_symbol == :r
                        path = CylindricalPoint.(path)
                        map(x -> x[2], path),
                        map(x -> x[3], path)
                    end
                    x, y
                end
            end
        end
    end
end
