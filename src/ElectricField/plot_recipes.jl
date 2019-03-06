function get_xyz_vector_from_rφz_field_vector_at_rφz(field,r,φ,z,ir,iφ,iz)::Vector
    startpoint_vector = get_xyz_vector_from_rφz_vector([r,φ,z])
    endpoint_vector = get_xyz_vector_from_rφz_vector([r,φ,z]+field[ir,iφ,iz])
    xyz_vector = endpoint_vector-startpoint_vector
    for ic in 1:size(xyz_vector,1)
        isapprox(xyz_vector[ic],0.0) ? xyz_vector[ic] = 0.0 : nothing
    end
    xyz_vector
end

function get_rφz_vector_from_xyz_vector(v)
    r = sqrt(v[1]^2 + v[2]^2)
    φ = atan(v[2],v[1])
    z = v[3]
    [r,φ,z]
end

function get_xy_magnitude(xyz_vector::AbstractArray)
    sqrt(xyz_vector[1]^2+xyz_vector[2]^2)
end

@recipe function f(electrical_field::Array{SVector{3,T},3}, grid::CylindricalGrid{T}; view=:Components, plane=:rz, i_fixed=3, spacing = 8, vectorscale = 0.0018, SI_factor=1/1000.) where T <: AbstractFloat

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
# @recipe function f(electrical_field::Array{SVector{3,T},3}, d::SolidStateDetector, grid::CylindricalGrid; φ_value=deg2rad(0), spacing=0.003, steps=1000, myscale=1, potential=true) where T <: AbstractFloat
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

@userplot myQuiver
@recipe function f(gdd::myQuiver; edges=0:1:3000)
    @series begin
        [1,2], [2,3]
    end
end


@recipe function f( electric_field::Array{ <:StaticVector{3, T}, 3}, det::SolidStateDetector, ep::ElectricPotential{T};
            φ=missing, spacing=0.003, n_steps=3000, potential=true, contours_equal_potential=true) where {T <: AbstractFloat}
    size --> (700,900)
    det.bulk_type == :ptype ? interpolated_efield = setup_interpolated_vectorfield(electric_field, ep.grid) : nothing
    det.bulk_type == :ntype ? interpolated_efield = setup_interpolated_vectorfield(map(x -> -1*x, electric_field), ep.grid) : nothing
    aspect_ratio --> 1
    title --> "Electric Field Lines @φ=$(φ)°"
    xlabel --> L"$r$ / m"
    ylabel --> L"$z$ / m"

    if ismissing(φ)
        φ = 0
    end
    φ_rad = deg2rad(φ)

    if potential==true
        @series begin
            contours_equal_potential --> contours_equal_potential
            φ --> φ
            ep
        end
    end

    corner_offset = 5e-5
    spawn_positions = []

    for (i,tuple) in enumerate(det.segmentation_r_ranges)
        if tuple[1]==tuple[2]
            for z in corner_offset + det.segmentation_z_ranges[i][1]:spacing:det.segmentation_z_ranges[i][2] - corner_offset
                push!(spawn_positions, MVector{3,T}(tuple[1],φ_rad,z))
            end
        end
    end

    for (i,tuple) in enumerate(det.segmentation_z_ranges)
        if tuple[1]==tuple[2]
            for r in corner_offset+det.segmentation_r_ranges[i][1]:spacing:det.segmentation_r_ranges[i][2] - corner_offset
                push!(spawn_positions, MVector{3,T}(r,φ_rad,tuple[1]))
            end
        end
    end

    for (i,orientation) in enumerate(det.segmentation_types)
        if orientation != "Tubs"
            if orientation[1]=='c' o=-1 else o=1 end
            for z in det.segmentation_z_ranges[i][1]:spacing:det.segmentation_z_ranges[i][2]
                push!(spawn_positions,MVector{3,T}(o * corner_offset+analytical_taper_r_from_φz(φ_rad, z, orientation, det.segmentation_r_ranges[i],
                                                                                                            det.segmentation_phi_ranges[i],
                                                                                                            det.segmentation_z_ranges[i]),
                                                                                                            φ_rad, z))
            end
        end
    end
    if typeof(det) <: Union{BEGe, Coax, InvertedCoax}
        for i in eachindex(spawn_positions)
            if spawn_positions[i][3] == geom_round(det.crystal_length)
                spawn_positions[i][3] -= 1e-5
            end
            if spawn_positions[i][1] == geom_round(det.crystal_radius)
                spawn_positions[i][1] -= 1e-5
            end
            if spawn_positions[i][3] == 0
                spawn_positions[i][3] += 1e-5
            end
            if typeof(det) <: Union{Coax, InvertedCoax}
                if spawn_positions[i][1] == geom_round(det.borehole_radius) #&& spawn_positions[i][3] != det.borehole_length
                    spawn_positions[i][1] += 1e-4
                end
                if spawn_positions[i][1] == geom_round(det.borehole_length )
                    spawn_positions[i][1] -= 1e-4
                end
            end
        end
    end
    spawn_positions_xyz::Vector{SVector{3, T}} = map(x -> CartesianPoint(CylindricalPoint(x[1],x[2],x[3])), spawn_positions)


    for i in eachindex(spawn_positions[1:end])
        path = [@SVector zeros(T,3) for i in 1:n_steps]
        drift_charge!(path, det, spawn_positions_xyz[i], T(1f-9), interpolated_efield)
        @series begin
            c --> :white
            label --> ""
            map(x->sqrt(x[1]^2+x[2]^2),path), map(x->x[3],path)
        end
    end
    # println(spawn_positions)
end
