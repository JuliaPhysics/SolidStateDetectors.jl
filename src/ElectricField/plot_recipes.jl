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
                grid[:r] ./ SI_factor, grid[:z] ./ SI_factor, vectorfield_r[:,i_fixed,:]'
            end
            @series begin
                subplot := 2
                colorbar_title --> "Field Strength [V / "*units[SI_factor]*"]"
                title := "φ_component"

                grid[:r] ./SI_factor, grid[:z] ./SI_factor, vectorfield_φ[:,i_fixed,:]'
            end
            @series begin
                subplot := 3
                title := "z_component"
                ylabel --> "z ["*units[SI_factor]*"]"
                xlabel --> "r ["*units[SI_factor]*"]"
                grid[:r] ./SI_factor, grid[:z] ./SI_factor, vectorfield_z[:,i_fixed,:]'
            end
            @series begin
                subplot := 4
                colorbar_title --> "Field Strength [V / "*units[SI_factor]*"]"
                xlabel --> "r ["*units[SI_factor]*"]"
                title:= "magnitude"
                grid[:r] ./SI_factor, grid[:z] ./SI_factor, vectorfield_magn[:,i_fixed,:]'
            end
        end
    elseif view == :ef
        if plane == :rφ
            vectorfield_xyz = Array{Vector{Float32}}(undef,size(vectorfield,1),size(vectorfield,2),size(vectorfield,3));
            for (iz,z) in enumerate(grid[:z])
                for (iφ,φ) in enumerate(grid[:φ])
                    for (ir,r) in enumerate(grid[:r])
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
        title := "z = $(round(grid[:z][i_fixed]/SI_factor,digits=2)) / mm"
        xlims := (-1.2/SI_factor*maximum(grid[:r]),1.2/SI_factor*maximum(grid[:r]))
        ylims := (-1.2/SI_factor*maximum(grid[:r]),1.2/SI_factor*maximum(grid[:r]))
        for (ir,r) in enumerate(grid[:r][1:spacing:end])
        for (iφ,φ) in enumerate(grid.φ)
                x= r*cos(φ)
            y= r*sin(φ)
            ir_actual=findfirst(x->x==r,grid[:r])
            iφ_actual=findfirst(x->x==φ,grid.φ)
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
#
# @recipe function f(electrical_field::Array{SVector{3,T},3}, d::SolidStateDetector, grid::CylindricalGrid; φ_value=deg2rad(0), spacing=0.003, steps=1000, myscale=1, potential=true) where{T <: SSDFloat}
#     size --> (900,1100)
#     d.bulk_type == :ptype ? interpolated_efield = setup_interpolated_vectorfield(electrical_field, grid) : nothing
#     d.bulk_type == :ntype ? interpolated_efield = setup_interpolated_vectorfield(map(x->-1*x,electrical_field), grid) : nothing
#     if potential==true
#         @series begin
#             φ --> φ_value
#             grid
#         end
#     end
#
#     corner_offset = 5e-5
#     spawn_positions::Array{Array{T,1},1}=[]
#
#     for (i,tuple) in enumerate(d.segmentation_r_ranges)
#         if tuple[1]==tuple[2]
#             for z in corner_offset + d.segmentation_z_ranges[i][1]:spacing:d.segmentation_z_ranges[i][2] - corner_offset
#                 push!(spawn_positions,[tuple[1],φ_value,z])
#             end
#         end
#     end
#
#     for (i,tuple) in enumerate(d.segmentation_z_ranges)
#         if tuple[1]==tuple[2]
#             for r in corner_offset+d.segmentation_r_ranges[i][1]:spacing:d.segmentation_r_ranges[i][2] - corner_offset
#                 push!(spawn_positions,[r,φ_value,tuple[1]])
#             end
#         end
#     end
#
#     for (i,orientation) in enumerate(d.segmentation_types)
#         if orientation != "Tubs"
#             if orientation[1]=='c' o=-1 else o=1 end
#             for z in d.segmentation_z_ranges[i][1]:spacing:d.segmentation_z_ranges[i][2]
#                 push!(spawn_positions,[o * corner_offset+analytical_taper_r_from_φz(φ_value, z, orientation, d.segmentation_r_ranges[i],
#                                                                                                             d.segmentation_phi_ranges[i],
#                                                                                                             d.segmentation_z_ranges[i]),
#                                                                                                             φ_value, z])
#             end
#         end
#     end
#
#     for i in eachindex(spawn_positions[1:end-1])
#
#         xpath, ypath, zpath = driftonecharge(d, get_xyz_vector_from_rφz_vector(spawn_positions[i]), steps, myscale*1e-9, :e, interpolated_efield, interpolated_efield)
#         rpath=[]
#
#         for ir in eachindex(xpath)
#             push!(rpath, sqrt(xpath[ir]^2+ypath[ir]^2))
#         end
#
#         @series begin
#             c --> :white
#             label --> ""
#             rpath, zpath
#         end
#     end
# end
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


# @recipe function f( electric_field::Array{ <:StaticVector{3, T}, 3}, det::SolidStateDetector, ep::ElectricPotential{T};
#             φ=missing, spacing=T(0.003), n_steps=3000, potential=true, contours_equal_potential=true, offset = T(5e-5)) where {T <: SSDFloat}
#     size --> (700,900)
#     det.bulk_type == :ptype ? interpolated_efield = setup_interpolated_vectorfield(electric_field, ep.grid) : nothing
#     det.bulk_type == :ntype ? interpolated_efield = setup_interpolated_vectorfield(map(x -> -1*x, electric_field), ep.grid) : nothing
#     if ismissing(φ)
#         φ = 0
#     end
#     φ_rad = deg2rad(φ)
#     aspect_ratio --> 1
#     title --> "Electric Field Lines @φ=$(φ)°"
#     xlabel --> L"$r$ / m"
#     ylabel --> L"$z$ / m"
#
#     if potential==true
#         @series begin
#             contours_equal_potential --> contours_equal_potential
#             φ --> φ
#             ep
#         end
#     end
#
#     spawn_positions_cyl::Vector{CylindricalPoint} = []
#
#     for contact in det.contacts
#         if (typeof(contact) == SSD.Contact{T,:N} && det.bulk_type == :ntype) || (typeof(contact) == SSD.Contact{T,:P} && det.bulk_type == :ptype)
#             nothing
#         else
#             for g in contact.geometry_positive
#                 if typeof(g) == SSD.Tube{T}
#                     rStart,rStop = get_r(g)
#                     zStart,zStop = get_z(g)
#                     for r in rStart:spacing:rStop
#                         push!(spawn_positions_cyl,CylindricalPoint{T}(r,φ,zStart+offset))
#                         push!(spawn_positions_cyl,CylindricalPoint{T}(r,φ,zStart-offset))
#                         push!(spawn_positions_cyl,CylindricalPoint{T}(r,φ,zStop+offset))
#                         push!(spawn_positions_cyl,CylindricalPoint{T}(r,φ,zStop-offset))
#                     end
#                     for z in zStart:spacing:zStop
#                         push!(spawn_positions_cyl,CylindricalPoint{T}(rStart+offset,φ,z))
#                         push!(spawn_positions_cyl,CylindricalPoint{T}(rStart-offset,φ,z))
#                         push!(spawn_positions_cyl,CylindricalPoint{T}(rStop+offset,φ,z))
#                         push!(spawn_positions_cyl,CylindricalPoint{T}(rStop-offset,φ,z))
#                     end
#                 elseif typeof(g) == SSD.ConeMantle{T}
#                     zStart,zStop = get_z(g)
#                     for z in zStart:spacing:zStop
#                         push!(spawn_positions_cyl,CylindricalPoint{T}(get_diagonal_r_from_z(g.cone, z)+offset,φ,z))
#                         push!(spawn_positions_cyl,CylindricalPoint{T}(get_diagonal_r_from_z(g.cone, z)-offset,φ,z))
#                         push!(spawn_positions_cyl,CylindricalPoint{T}(get_diagonal_r_from_z(g.cone, z)+offset,φ,z))
#                         push!(spawn_positions_cyl,CylindricalPoint{T}(get_diagonal_r_from_z(g.cone, z)-offset,φ,z))
#                     end
#                 end
#             end
#         end
#     end
#     filter!(x -> x in det && !in(x, det.contacts),spawn_positions_cyl)
#     spawn_positions_xyz::Vector{CartesianPoint{T}} = map(x -> CartesianPoint(CylindricalPoint{T}(x[1],x[2],x[3])), spawn_positions_cyl)
#     for i in eachindex(spawn_positions_cyl[1:end])
#         path = [@SVector zeros(T,3) for i in 1:n_steps]
#         drift_charge!(path, det, spawn_positions_xyz[i], T(1f-9), interpolated_efield)
#         @series begin
#             c --> :white
#             label --> ""
#             map(x->sqrt(x[1]^2+x[2]^2),path), map(x->x[3],path)
#         end
#     end
# end

@userplot Plot_electric_field
@recipe function f(gdd::Plot_electric_field; φ = missing, spacing = 10, grid_spacing=[0.0005, deg2rad(1.0), 0.0005], n_steps=3000, potential=true, contours_equal_potential=true, offset = (5e-5))
    setup = gdd.args[1]
    T = Float32
    φ = ismissing(φ) ? T(0) : T(φ)
    φ_rad = deg2rad(φ)
    aspect_ratio --> 1
    title --> "Electric Field Lines @φ=$(round(φ, digits=2))°"
    xlabel --> L"$r$ / m"
    ylabel --> L"$z$ / m"

    if potential==true
        @series begin
            contours_equal_potential --> contours_equal_potential
            φ --> φ
            setup.electric_potential
        end
    end


    contacts_to_spawn_charges_for = filter!(x -> x.id !=1, SSD.Contact{T}[c for c in setup.detector.contacts])
    spawn_positions = CylindricalPoint{T}[]
    grid = Grid(setup.detector, init_grid_spacing = grid_spacing, full_2π = true)
    pt_offset = T[offset,0.0,offset]

    for c in contacts_to_spawn_charges_for
        ongrid_positions= map(x-> CylindricalPoint{T}(grid[x...]),SSD.paint_object(c, grid, φ_rad))
        for position in ongrid_positions
            push!(spawn_positions, CylindricalPoint{T}((position + pt_offset)...))
            push!(spawn_positions, CylindricalPoint{T}((position - pt_offset)...))
        end
    end

    filter!(x -> x in setup.detector && !in(x, setup.detector.contacts), spawn_positions)

    el_field_itp = SSD.get_interpolated_drift_field(setup.electric_field.data, setup.electric_field.grid)
    el_field_itp_inv = SSD.get_interpolated_drift_field(setup.electric_field.data .* -1, setup.electric_field.grid)
    for (ipos, pos) in enumerate(spawn_positions)
        if ((spacing-1)+ipos)%spacing == 0
            path = CartesianPoint{T}[CartesianPoint{T}(0.0,0.0,0.0) for i in 1:n_steps]
            SSD.drift_charge!(path, setup.detector, setup.point_types, setup.electric_potential.grid, CartesianPoint(pos), T(2e-9), el_field_itp, verbose = false )
            @series begin
                c --> :white
                label --> ""
                map(x->sqrt(x[1]^2+x[2]^2),path), map(x->x[3],path)
            end

            path = CartesianPoint{T}[CartesianPoint{T}(0.0,0.0,0.0) for i in 1:n_steps]
            SSD.drift_charge!(path, setup.detector, setup.point_types, setup.electric_potential.grid, CartesianPoint(pos), T(2e-9), el_field_itp_inv, verbose = false )
            @series begin
                c --> :white
                label --> ""
                map(x->sqrt(x[1]^2+x[2]^2),path), map(x->x[3],path)
            end

        end
    end
end
