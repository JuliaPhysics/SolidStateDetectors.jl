@recipe function f(det::SolidStateDetector{T}; n = 30, φ = missing, seriescolor = missing, label = missing) where {T}

    clabel = (ismissing(label) ? map(c -> (c.name == "" ? c.id : c.name), det.contacts) : label)
    if !(typeof(clabel) <: AbstractArray) clabel = [clabel] end
    ccolor = (ismissing(seriescolor) ? map(c -> c.id, det.contacts) : seriescolor)
    if !(typeof(ccolor) <: AbstractArray) ccolor = [ccolor] end

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


@recipe function f(contact::Contact{T}; n = 30, seriescolor = missing) where {T}
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
            c
        end
    end
end
