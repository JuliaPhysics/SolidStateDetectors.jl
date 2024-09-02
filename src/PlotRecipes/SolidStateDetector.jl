@recipe function f(p::Passive{T}) where {T}
    seriestype --> :csg
    linecolor --> :silver
    seriescolor --> :silver
    st = plotattributes[:seriestype]
    if st in [:x, :y, :z]
        n_samples --> 200
    end
    fillalpha --> 0.2
    l = p.name != "" ? p.name : "Passive $(p.id)"
    @series begin
        #add empty line so that label is shown with alpha = 1
        seriestype := :path
        label --> l
        linewidth := 1
        if st in [:x, :y, :z]
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
    if st in [:x, :y, :z]
        n_samples --> 200
    end
    @series begin
        #add empty line so that label is shown with alpha = 1
        seriestype := :path
        label --> "Semiconductor"
        linewidth := 1
        if st in [:x, :y, :z]
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
    if st in [:x, :y, :z]
        n_samples --> 200
    end
    fillalpha --> 0.2
    l = contact.name != "" ? "$(contact.name) (id: $(contact.id))" : "Contact - id: $(contact.id)"
    @series begin
        #add empty line so that label is shown with alpha = 1
        seriestype := :path
        label --> l
        linewidth := 1
        if st in [:x, :y, :z]
            T[]*internal_length_unit, T[]*internal_length_unit
        else
            T[]*internal_length_unit, T[]*internal_length_unit, T[]*internal_length_unit
        end
    end 
    label := ""
    contact.geometry
end

@recipe function f(det::SolidStateDetector; show_semiconductor = false, show_passives = true, n_samples = 40)
    seriestype --> :csg
    show_normal --> false
    st = plotattributes[:seriestype]
    if st in [:x, :y, :z]
        n_samples --> 200
    end
    
    plot_objects = []
    append!(plot_objects, det.contacts)
    if show_semiconductor
        push!(plot_objects, det.semiconductor)
    end
    if show_passives && !ismissing(det.passives)
        append!(plot_objects, det.passives)
    end
    
    if haskey(plotattributes, :seriestype) && plotattributes[:seriestype] == :samplesurface
        scales = [get_scale(o.geometry) for o in plot_objects]
        max_scale = maximum(scales)
    end
    
    for (i,o) in enumerate(plot_objects)
        @series begin
            if haskey(plotattributes, :seriestype) && plotattributes[:seriestype] == :samplesurface
                n_samples --> Int(ceil(n_samples * scales[i]/max_scale))
                CSG_scale := scales[i]
            end
            o
        end
    end
end
