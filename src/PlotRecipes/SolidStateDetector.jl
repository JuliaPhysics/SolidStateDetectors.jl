@recipe function f(det::SolidStateDetector{T}; SSD_style = :wireframe, n = 30, φ = missing, seriescolor = missing, label = missing, alpha_factor = 1) where {T}
    if !(SSD_style in [:wireframe])#, :samplesurface])
        #@warn "Chose SSD_style from [:wireframe, :samplesurface]. Defaulting to :wireframe"
        SSD_style = :wireframe
    end
    clabel = (ismissing(label) ? map(c -> (c.name == "" ? c.id : c.name), det.contacts) : label)
    if !(typeof(clabel) <: AbstractArray) clabel = [clabel] end
    ccolor = (ismissing(seriescolor) ? map(c -> c.id, det.contacts) : seriescolor)
    if !(typeof(ccolor) <: AbstractArray) ccolor = [ccolor] end
    world_size = missing
    #=
    if SSD_style == :samplesurface
        grid = Grid(det) # DOES NOT WORK ANYMORE
        CS = get_coordinate_system(det) # DOES NOT WORK ANYMORE
        if CS == Cylindrical
            world_size = CylindricalVector{T}(width(grid.r.interval), π, width(grid.z.interval))
        elseif CS == Cartesian
            world_size = CartesianVector{T}(width(grid.x.interval), width(grid.y.interval), width(grid.z.interval))
        else
            @error "Could not determine the world size, try SSD_style = :wireframe"
        end
    end
    =#
    if ismissing(φ)
        xguide --> "x / m"
        yguide --> "y / m"
        zguide --> "z / m"

        for (cn,contact) in enumerate(det.contacts)
            linewidth --> 2
            seriescolor := ccolor[(cn-1)%length(ccolor)+1]
            @series begin
                label := ""
                n --> n
                SSD_style --> SSD_style
                world_size --> world_size
                alpha_factor --> alpha_factor
                contact
            end
            @series begin
                label := clabel[(cn-1)%length(clabel)+1]
                seriescolor := ccolor[(cn-1)%length(ccolor)+1]
                []
            end
        end
    else
        @info "2D plotting of the detector not yet implemented."
    end
end


@recipe function f(contact::Contact{T}; SSD_style = :wireframe, n = 30, seriescolor = missing, world_size = missing, alpha_factor = 1) where {T}
    ccolor = (ismissing(seriescolor) ? contact.id : seriescolor)
    @series begin
        seriescolor := ccolor
        label --> (contact.name == "" ? contact.id : contact.name)
        []
    end
    for (cn,c) in enumerate(contact.geometry_positive)
        @series begin
            seriescolor := ccolor
            label := ""
            n --> n
            SSD_style --> SSD_style
            world_size --> world_size
            geometry_negative --> contact.geometry_negative
            alpha_factor --> alpha_factor
            c
        end
    end
end