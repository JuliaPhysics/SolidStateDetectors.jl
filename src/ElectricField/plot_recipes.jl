function get_xyz_vector_from_rφz_field_vector_at_rφz(field,r,φ,z,ir,iφ,iz)::Vector
    startpoint_vector = get_xyz_vector_from_rφz_vector([r,φ,z])
    endpoint_vector = get_xyz_vector_from_rφz_vector([r,φ,z]+field[ir,iφ,iz])
    xyz_vector = endpoint_vector-startpoint_vector
    for ic in 1:size(xyz_vector,1)
        isapprox(xyz_vector[ic],0.0) ? xyz_vector[ic] = 0.0 : nothing
    end
    xyz_vector
end


function get_xy_magnitude(xyz_vector::AbstractArray)
    sqrt(xyz_vector[1]^2+xyz_vector[2]^2)
end

@recipe function f(electrical_field::Array{SVector{3,T},3}, grid::CylindricalGrid{T}; view=:Components, plane=:rz, i_fixed=3, spacing = 8, vectorscale = 0.0018, SI_factor=1/1000.) where{T <: SSDFloat}

    vectorfield = electrical_field.*SI_factor
    units = Dict(1e-3=>"mm",1e-2=>"cm",1e-1=>"dm",1=>"m",1.0=>"m")
    if view == :Components
        if plane == :rz
            layout := (2,2)
            size := (800,900)
            vectorfield_r = get_component_field(vectorfield);
            vectorfield_φ = get_component_field(vectorfield,:phi);
            vectorfield_z = get_component_field(vectorfield,:z);
            vectorfield_magn = get_magnitude_of_rφz_vector.(vectorfield);

            st := :heatmap
            # colorbar_title := "Field Strength / V m\$\^\{-1\}\$"
            @series begin
                subplot := 1
                title := "r_component"
                ylabel --> "z ["*units[SI_factor]*"]"
                grid.axes[1] ./ SI_factor, grid.axes[3] ./ SI_factor, vectorfield_r[:,i_fixed,:]'
            end
            @series begin
                subplot := 2
                colorbar_title --> "Field Strength [V / "*units[SI_factor]*"]"
                title := "φ_component"

                grid.axes[1] ./SI_factor, grid.axes[3] ./SI_factor, vectorfield_φ[:,i_fixed,:]'
            end
            @series begin
                subplot := 3
                title := "z_component"
                ylabel --> "z ["*units[SI_factor]*"]"
                xlabel --> "r ["*units[SI_factor]*"]"
                grid.axes[1] ./SI_factor, grid.axes[3] ./SI_factor, vectorfield_z[:,i_fixed,:]'
            end
            @series begin
                subplot := 4
                colorbar_title --> "Field Strength [V / "*units[SI_factor]*"]"
                xlabel --> "r ["*units[SI_factor]*"]"
                title:= "magnitude"
                grid.axes[1] ./SI_factor, grid.axes[3] ./SI_factor, vectorfield_magn[:,i_fixed,:]'
            end
        end
    elseif view == :ef
        if plane == :rφ
            vectorfield_xyz = Array{Vector{T}}(undef,size(vectorfield,1),size(vectorfield,2),size(vectorfield,3));
            for (iz,z) in enumerate(grid.axes[3])
                for (iφ,φ) in enumerate(grid.axes[2])
                    for (ir,r) in enumerate(grid.axes[1])
                        vectorfield_xyz[ir,iφ,iz]=get_xyz_vector_from_rφz_field_vector_at_rφz(vectorfield,r,φ,z,ir,iφ,iz)
                    end
                end
            end
        elseif plane == :xy
            vectorfield_xyz = electrical_field.field
        end
        vectorfield_xy_magn = map(x->get_xy_magnitude(x),vectorfield_xyz[:,:,i_fixed])
        max_magn = maximum(vectorfield_xy_magn)
        diff_magn = max_magn-minimum(vectorfield_xy_magn)

        size := (800,800)
        line := (:arrow,:blue,2)
        label := ""
        ylabel := "y "
        xlabel := "x "
        title := "z = $(round(grid.axes[3][i_fixed]/SI_factor,digits=2)) / mm"
        xlims := (-1.2/SI_factor*maximum(grid.axes[1]),1.2/SI_factor*maximum(grid.axes[1]))
        ylims := (-1.2/SI_factor*maximum(grid.axes[1]),1.2/SI_factor*maximum(grid.axes[1]))
        for (ir,r) in enumerate(grid.axes[1][1:spacing:end])
        for (iφ,φ) in enumerate(grid..axes[2])
                x= r*cos(φ)
            y= r*sin(φ)
            ir_actual=findfirst(x->x==r,grid.axes[1])
            iφ_actual=findfirst(x->x==φ,grid.axes[2])
                xy_magn = vectorfield_xy_magn[ir_actual,iφ_actual]
                vector=vectorfield_xyz[ir_actual,iφ_actual,i_fixed]/xy_magn
                vector*=((vectorfield_xy_magn[ir_actual,iφ_actual]-0.8*minimum(vectorfield_xy_magn))/diff_magn)
                vector*=vectorscale
                @series begin
                    [x-0.5*vector[1],x+0.5*vector[1]]/SI_factor, [y-0.5*vector[2],y+0.5*vector[2]]/SI_factor
                end
            end
        end
    end
end

using Base.Math
@recipe function f( electric_field::Array{ <:StaticVector{3, T}, 3}, det::SolidStateDetector; φ=missing, z = missing, scaling_factor=1.0, spacings = [3,20]) where T
    ismissing(φ) && !ismissing(z) ? iz = searchsortednearest(det.zs,T(z)) : nothing
    !ismissing(φ) && ismissing(z) ? iφ = searchsortednearest(det.φs,T(φ)) : nothing
    guidefontsize --> 15
    tickfontsize -->15
    # labelfontsize --> 13
    if ismissing(φ)
        xy = T[0, 0]
        values = T[0, 0]
        magnitudes_xy = []
        for (iφ ,φ)  in enumerate(det.φs)
            if (iφ+spacings[1]-1)%spacings[1]==0
                for (ir, r) in enumerate(det.rs)
                    if (ir+spacings[2]-1)%spacings[2] == 0
                        xy = hcat(xy,[ r * cos(φ), r * sin(φ)])
                        values = hcat(values, electric_field[ir, iφ, iz][1:2])
                        push!(magnitudes_xy,sqrt(electric_field[ir, iφ, iz][1]^2 + electric_field[ir, iφ, iz][2]^2))
                    end
                end
            end
        end

        num = minimum([mean(diff(det.φs)*spacings[1]),mean(diff(det.rs))*spacings[2]])
        nom = maximum(magnitudes_xy)
        scaling = num/nom * scaling_factor
        for i in 1:size(xy,2)
            @series begin
                title --> "z = $(round(det.zs[iz] .* 1000,digits=2)) mm"
                xlabel --> "x / mm"
                ylabel --> "y / mm"
                aspect_ratio --> 1
                label --> ""
                color --> :red
                arrow --> true
                [xy[1,i] - 0.5*scaling*values[1,i], xy[1,i] + 0.5*scaling*values[1,i]] .* 1000, [xy[2,i] - 0.5*scaling*values[2,i], xy[2,i] + 0.5*scaling*values[2,i]] .* 1000
            end
        end

    elseif ismissing(z)
        rz = T[0, 0]
        values = T[0, 0]
        magnitudes_rz = []
        for (iz ,z)  in enumerate(det.zs)
            if (iz+spacings[1])%spacings[1]==0
                for (ir, r) in enumerate(det.rs)
                    if (ir+spacings[2]-1)%spacings[2] == 0
                        rz = hcat(rz,[ r, z])
                        values = hcat(values, electric_field[ir, iφ, iz][1:2:3])
                        push!(magnitudes_rz,sqrt(electric_field[ir, iφ, iz][1]^2 + electric_field[ir, iφ, iz][3]^2))
                    end
                end
            end
        end
        num = minimum([mean(diff(det.zs)*spacings[1]),mean(diff(det.rs))*spacings[2]])
        nom = maximum(magnitudes_rz)
        scaling = num/nom * scaling_factor
        for i in 1:size(rz,2)
            @series begin
                title --> "φ = $(round(rad2deg(det.φs[iφ]) , digits=2)) deg"
                xlabel --> "r [mm]"
                ylabel --> "z [mm]"
                aspect_ratio --> 1
                label --> ""
                color --> :red
                arrow --> true
                [rz[1,i] - 0.5*scaling*values[1,i], rz[1,i] + 0.5*scaling*values[1,i]] .* 1000, [rz[2,i] - 0.5*scaling*values[2,i], rz[2,i] + 0.5*scaling*values[2,i]] .* 1000
            end
        end
    end
end

@userplot MyQuiver
@recipe function f(gdd::MyQuiver; scaling=1.0)
    xy::Matrix = gdd.args[1]
    values = gdd.args[2]
    for i in 1:size(xy,2)
        @series begin
            arrow --> true
            [xy[1,i] - 0.5*scaling*values[1,i], xy[1,i] + 0.5*scaling*values[1,i]], [xy[2,i] - 0.5*scaling*values[2,i], xy[2,i] + 0.5*scaling*values[2,i]]
        end
    end
end

@recipe function f( ef::ElectricField{T, 3, :cylindrical};
                    r = missing,
                    φ = missing,
                    z = missing,
                    contours_equal_potential=false,
                    full_det = false ) where {T}
    g::Grid{T, 3, :cylindrical} = ef.grid

    seriescolor --> :inferno
    st --> :heatmap
    aspect_ratio --> 1
    foreground_color_border --> nothing
    tick_direction --> :out

    ef_magn = get_magnitude_of_rφz_vector.(ef)

    cross_section::Symbol, idx::Int, idx_mirror = if ismissing(φ) && ismissing(r) && ismissing(z)
        :φ, 1, Int(length(g.φ)/2)
    elseif !ismissing(φ) && ismissing(r) && ismissing(z)
        φ_rad::T = T(deg2rad(φ))
        while !(g.φ.interval.left <= φ_rad <= g.φ.interval.right) && g.φ.interval.right != g.φ.interval.left
            if φ_rad > g.φ.interval.right
                φ_rad -= g.φ.interval.right - g.φ.interval.left
            elseif φ_rad < g.φ.interval.left
                φ_rad += g.φ.interval.right - g.φ.interval.left
            end
        end
        :φ, searchsortednearest(g.φ, φ_rad), searchsortednearest(g.φ, (φ_rad+π)%(2π))
    elseif ismissing(φ) && !ismissing(r) && ismissing(z)
        :r, searchsortednearest(g.r, T(r)), searchsortednearest(g.r, T(r))
    elseif ismissing(φ) && ismissing(r) && !ismissing(z)
        :z, searchsortednearest(g.z, T(z)), searchsortednearest(g.z, T(z))
    else
        error(ArgumentError, ": Only one of the keywords `r, φ, z` is allowed.")
    end
    value::T = if cross_section == :φ
        g.φ[idx]
    elseif cross_section == :r
        g.r[idx]
    elseif cross_section == :z
        g.z[idx]
    end

    @series begin
        if cross_section == :φ
            title --> "Electric Potential @$(cross_section) = $(round(rad2deg(value), sigdigits = 2))"
            xlabel --> "r / m"
            ylabel --> "z / m"
            if full_det == true
                size --> ( 400, 350 / (g.r[end] - g.r[1]) * (g.z[end] - g.z[1]) )
                vcat(-1 .* g.r[end:-1:2], g.r),  g.z, cat(ef_magn[end:-1:2, idx_mirror, :]', ef_magn[:, idx, :]', dims = 2)
            else
                size --> ( 400, 350 / (g.r[end] - g.r[1]) * (g.z[end] - g.z[1]) )
                g.r, g.z, ef_magn[:, idx,:]'
            end
        elseif cross_section == :r
            title --> "Electric Potential @$(cross_section) = $(round(value, sigdigits = 2))"
            g.φ, g.z, ef_magn[idx,:,:]'
        elseif cross_section == :z
            title --> "Electric Potential @$(cross_section) = $(round(value, sigdigits = 2))"
            proj --> :polar
            g.φ, g.r, ef_magn[:,:,idx]
        end
    end
    if contours_equal_potential
        @series begin
            seriescolor := :thermal
            st := :contours
            if cross_section == :φ
                g.r, g.z, ef_magn[:, idx,:]'
            elseif cross_section == :r
                g.φ, g.z, ef_magn[idx,:,:]'
            elseif cross_section == :z
                proj --> :polar
                g.φ, g.r, ef_magn[:,:,idx]
            end
        end
    end
end


@userplot Plot_electric_field_
@recipe function f(gdd::Plot_electric_field_; φ = missing, r = missing, x = missing, y = missing, z = missing,
                    potential = false, contours_equal_potential = true, full_det = false, field_lines = true,
                    sampling = 4u"mm", spacing = 1, max_nsteps = 3000, offset = (0.001)
                    )
    sim = gdd.args[1]
    S = get_coordinate_system(sim.electric_field.grid)
    T = SolidStateDetectors.get_precision_type(sim.detector)

    dim_array = [φ, r, x, y, z]
    dim_symbols_array = [:φ, :r, :x, :y, :z]
    if isempty(skipmissing(dim_array))
        if S == :cylindrical
            v::T = 0
            dim_number = 2
            dim_symbol = :φ
        elseif S == :cartesian
            dim_number = 1
            dim_symbol = :x
            v = sim.electric_field.grid[dim_number][ div(length(sim.electric_field.grid[dim_number]), 2) ]
        end
    elseif sum(ismissing.(dim_array)) == 4
        dim_number = findfirst(x -> !ismissing(x), dim_array)
        dim_symbol = dim_symbols_array[dim_number]
        v = dim_array[dim_number]
        if dim_symbol == :φ v = deg2rad(v) end
    else
        throw(ArgumentError("Only one keyword for a certain dimension is allowed. One of 'φ', 'r', 'x', 'y', 'z'"))
    end

    S::Symbol = get_coordinate_system(sim.detector)
    aspect_ratio --> 1
    title --> (field_lines ? "Electric Field Lines @$(dim_symbol)=$(round(dim_symbol == :φ ? rad2deg(v) : v, digits=2))" : "Electric Field @$(dim_symbol)=$(round(dim_symbol == :φ ? rad2deg(v) : v, digits=2))")
    xlabel --> (S == :cylindrical ? (dim_symbol == :r ? L"$\varphi$ / °" : L"$r$ / m") : (dim_symbol == :x ? L"$y$ / m" : L"$x$ / m"))
    ylabel --> L"$z$ / m"
        @series begin
            # contours_equal_potential --> contours_equal_potential
            if dim_symbol == :φ
                φ --> v
                full_det --> full_det
            elseif dim_symbol == :z
                if S == :cylindrical error("Electric field plot not yet implemented for z in cylindrical coordinates") end
                z --> v
                if S == :cylindrical
                    proj --> :polar
                end
            elseif dim_symbol == :r
                error("Electric field plot not yet implemented for r")
                r --> v
            elseif dim_symbol == :x
                x --> v
            elseif dim_symbol == :y
                y --> v
            end
            #clims --> missing
            potential ? sim.electric_potential : sim.electric_field
        end
end

@userplot Plot_electric_fieldlines
@recipe function f(gdd::Plot_electric_fieldlines; φ = missing, r = missing, x = missing, y = missing, z = missing,
                    sampling = 4u"mm", spacing = 2, max_nsteps=3000,
                    potential=true, contours_equal_potential=true, offset = (0.002), full_det = false)
    sim = gdd.args[1]
    S = get_coordinate_system(sim.electric_field.grid)
    T = SolidStateDetectors.get_precision_type(sim.detector)

    dim_array = [φ, r, x, y, z]
    dim_symbols_array = [:φ, :r, :x, :y, :z]
    if isempty(skipmissing(dim_array))
        if S == :cylindrical
            v::T = 0
            dim_number = 2
            dim_symbol = :φ
        elseif S == :cartesian
            dim_number = 1
            dim_symbol = :x
            v = sim.electric_field.grid[dim_number][ div(length(sim.electric_field.grid[dim_number]), 2) ]
        end
    elseif sum(ismissing.(dim_array)) == 4
        dim_number = findfirst(x -> !ismissing(x), dim_array)
        dim_symbol = dim_symbols_array[dim_number]
        v = dim_array[dim_number]
        if dim_symbol == :φ v = deg2rad(v) end
    else
        throw(ArgumentError("Only one keyword for a certain dimension is allowed. One of 'φ', 'r', 'x', 'y', 'z'"))
    end

    S::Symbol = get_coordinate_system(sim.detector)
    aspect_ratio --> 1
    title --> "Electric Field Lines @$(dim_symbol)=$(round(dim_symbol == :φ ? rad2deg(v) : v, digits=2))"
    xlabel --> (S == :cylindrical ? (dim_symbol == :r ? L"$\varphi$ / °" : L"$r$ / m") : (dim_symbol == :x ? L"$y$ / m" : L"$x$ / m"))
    ylabel --> L"$z$ / m"

    PT = S == :cylindrical ? CylindricalPoint{T} : CartesianPoint{T}

    contacts_to_spawn_charges_for = filter!(x -> x.id !=1, Contact{T}[c for c in sim.detector.contacts])
    spawn_positions = CartesianPoint{T}[]
    spawn_positions_mirror = CartesianPoint{T}[]
    grid = sim.electric_field.grid
    pt_offset = CartesianVector(CartesianPoint(CylindricalPoint{T}(offset, φ, offset)))
    pt_offset_mirror = CartesianVector(CartesianPoint(CylindricalPoint{T}(offset, φ+π, offset)))
    sampling_vector = T.(ustrip.([uconvert(u"m", sampling) for i in 1:3]))
    sampling_vector[2] = 2π
    for c in contacts_to_spawn_charges_for[:]#[4:4]
        for positive_geometry in c.geometry_positive[:]#[6:6]
            sampled_point_cyl = CylindricalPoint{T}.(SolidStateDetectors.sample(positive_geometry, sampling_vector))
            sampled_point_cyl =  map(x->geom_round(CylindricalPoint{T}(x[1], φ, x[3])), sampled_point_cyl)
            if sampled_point_cyl[1] in positive_geometry
                push!(spawn_positions, broadcast(+, CartesianPoint.(sampled_point_cyl), pt_offset )...)
                push!(spawn_positions, broadcast(-, CartesianPoint.(sampled_point_cyl), pt_offset )...)
            end
            if full_det == true
                sampled_point_cyl_mirror = map(x->geom_round(CylindricalPoint{T}(x[1], φ+π, x[3])), sampled_point_cyl)
                if sampled_point_cyl_mirror[1] in positive_geometry
                    push!(spawn_positions_mirror, broadcast(+, CartesianPoint.(sampled_point_cyl_mirror), pt_offset_mirror )...)
                    push!(spawn_positions_mirror, broadcast(-, CartesianPoint.(sampled_point_cyl_mirror), pt_offset_mirror )...)
                end
            end
        end
    end


    filter!(x -> x in sim.detector && !in(x, sim.detector.contacts), spawn_positions)
    full_det == true ? filter!(x ->x in sim.detector && !in(x, sim.detector.contacts), spawn_positions_mirror) : nothing
    spawn_positions = vcat(spawn_positions,spawn_positions_mirror)
    @info "$(round(Int,length(spawn_positions)/spacing)) drifts are now being simulated..."

    el_field_itp     = get_interpolated_drift_field(sim.electric_field.data,       sim.electric_field.grid)
    el_field_itp_inv = get_interpolated_drift_field(sim.electric_field.data .* -1, sim.electric_field.grid)

    @showprogress for (ipos, pos) in enumerate(spawn_positions)
        if ((spacing-1)+ipos)%spacing == 0

            path = CartesianPoint{T}[CartesianPoint{T}(0.0,0.0,0.0) for i in 1:max_nsteps]
            _drift_charge!(path, Vector{T}(undef, max_nsteps), sim.detector, sim.point_types, sim.electric_potential.grid, CartesianPoint(pos), T(2e-9), el_field_itp, verbose = false )
            filter!(x->x != CartesianPoint{T}(0.0,0.0,0.0), path)
            @series begin
                c --> :white
                if dim_symbol == :z && S == :cylindrical proj --> :polar end
                label --> ""
                x, y = if dim_symbol == :φ
                    map(x -> (pos in spawn_positions_mirror ? -1 : 1) * sqrt(x[1]^2+x[2]^2), path),
                    map(x -> x[3], path)
                elseif dim_symbol == :x
                    map(x -> x[2], path),
                    map(x -> x[3], path)
                elseif dim_symbol == :y
                    map(x -> x[1], path),
                    map(x -> x[3], path)
                elseif dim_symbol == :z
                    map(x -> x[1], path),
                    map(x -> x[2], path)
                end
                x, y
            end

            path = CartesianPoint{T}[CartesianPoint{T}(0.0,0.0,0.0) for i in 1:max_nsteps]
            _drift_charge!(path, Vector{T}(undef, max_nsteps), sim.detector, sim.point_types, sim.electric_potential.grid, CartesianPoint(pos), T(2e-9), el_field_itp_inv, verbose = false )
            filter!(x->x != CartesianPoint{T}(0.0,0.0,0.0), path)
            @series begin
                c --> :white
                if dim_symbol == :z && S == :cylindrical proj --> :polar end
                label --> ""
                x, y = if dim_symbol == :φ
                    map(x -> (pos in spawn_positions_mirror ? -1 : 1) * sqrt(x[1]^2+x[2]^2), path),
                    map(x -> x[3], path)
                elseif dim_symbol == :x
                    map(x -> x[2], path),
                    map(x -> x[3], path)
                elseif dim_symbol == :y
                    map(x -> x[1], path),
                    map(x -> x[3], path)
                elseif dim_symbol == :z
                    map(x -> x[1], path),
                    map(x -> x[2], path)
                end
                x, y
            end

        end
    end
end
export plot_electric_field
function plot_electric_field(sim; field_lines = true, kwargs...)
    plot_electric_field_(sim; kwargs...)
    if field_lines == true
        plot_electric_fieldlines!(sim; kwargs...)
    end
end
