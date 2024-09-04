@recipe function f(p::Passive{T}) where {T}
    seriestype --> :csg
    linecolor --> :silver
    seriescolor --> :silver
    st = plotattributes[:seriestype]
    fillalpha --> 0.2
    l = p.name != "" ? p.name : "Passive $(p.id)"
    @series begin
        #add empty line so that label is shown with alpha = 1
        seriestype := :path
        label --> l
        linewidth := 1
        if st == :slice
            T[]*internal_length_unit, T[]*internal_length_unit
        else
            T[]*internal_length_unit, T[]*internal_length_unit, T[]*internal_length_unit
        end
    end 
    label := ""
    p.geometry
end

@recipe function f(sc::Semiconductor{T}) where {T}
    seriestype --> :csg
    linecolor --> :grey
    seriescolor --> :grey
    st = plotattributes[:seriestype]
    @series begin
        #add empty line so that label is shown with alpha = 1
        seriestype := :path
        label --> "Semiconductor"
        linewidth := 1
        if st == :slice
            T[]*internal_length_unit, T[]*internal_length_unit
        else
            T[]*internal_length_unit, T[]*internal_length_unit, T[]*internal_length_unit
        end
    end 
    label := ""
    sc.geometry
end

@recipe function f(contact::Contact{T}) where {T}
    seriestype --> :csg
    seriescolor --> contact.id
    linecolor --> contact.id
    st = plotattributes[:seriestype]
    fillalpha --> 0.2
    l = contact.name != "" ? "$(contact.name) (id: $(contact.id))" : "Contact - id: $(contact.id)"
    @series begin
        #add empty line so that label is shown with alpha = 1
        seriestype := :path
        label --> l
        linewidth := 1
        if st == :slice
            T[]*internal_length_unit, T[]*internal_length_unit
        else
            T[]*internal_length_unit, T[]*internal_length_unit, T[]*internal_length_unit
        end
    end 
    label := ""
    contact.geometry
end

@recipe function f(det::SolidStateDetector; show_semiconductor = false, show_passives = true)
    seriestype --> :csg
    show_normal --> false
    st = plotattributes[:seriestype]

    plot_objects = []
    append!(plot_objects, det.contacts)
    if show_semiconductor
        push!(plot_objects, det.semiconductor)
    end
    if show_passives && !ismissing(det.passives)
        append!(plot_objects, det.passives)
    end
    
    max_scale = st in (:samplesurface, :slice) ? maximum(broadcast(o -> get_scale(o.geometry), plot_objects)) : missing
    
    for o in plot_objects
        @series begin
            CSG_scale := max_scale
            o
        end
    end
end
