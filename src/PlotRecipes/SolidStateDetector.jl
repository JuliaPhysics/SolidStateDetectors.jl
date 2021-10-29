@recipe function f(p::Passive{T}) where {T}
    seriestype --> :csg
    linecolor --> :silver
    #fillcolor --> :grey
    seriescolor --> :silver
    fillalpha --> 0.2
    l = p.name != "" ? p.name : "Passive $(p.id)"
    @series begin
        #add empty line so that label is shown with alpha = 1
        seriestype := :path
        label --> l
        linewidth := 1
        T[], T[]
    end 
    label := ""
    p.geometry
end

@recipe function f(sc::Semiconductor{T}) where {T}
    seriestype --> :csg
    linecolor --> :grey
    seriescolor --> :grey
    @series begin
        #add empty line so that label is shown with alpha = 1
        seriestype := :path
        label --> "Semiconductor"
        linewidth := 1
        T[], T[]
    end 
    label := ""
    sc.geometry
end

@recipe function f(contact::Contact{T}) where {T}
    seriestype --> :csg
    seriescolor --> contact.id
    linecolor --> contact.id
    fillcolor --> contact.id
    fillalpha --> 0.2
    l = contact.name != "" ? "$(contact.name) (id: $(contact.id))" : "Contact - id: $(contact.id)"
    @series begin
        #add empty line so that label is shown with alpha = 1
        seriestype := :path
        label --> l
        linewidth := 1
        T[], T[]
    end 
    label := ""
    contact.geometry
end

@recipe function f(det::SolidStateDetector; show_semiconductor = false, show_passives = true, n_samples = 100)

    show_normal --> false
    
    plot_objects = []
    append!(plot_objects, det.contacts)
    if show_semiconductor
        push!(plot_objects, det.semiconductor)
    end
    if show_passives && !ismissing(det.passives)
        append!(plot_objects, det.passives)
    end
    
    scales = nothing
    max_scale = nothing
    if haskey(plotattributes, :seriestype) && plotattributes[:seriestype] == :samplesurface
        scales = [get_scale(o.geometry) for o in plot_objects]
        max_scale = maximum(scales)
    end
    
    for (i,o) in enumerate(plot_objects)
        @series begin
            if haskey(plotattributes, :seriestype) && plotattributes[:seriestype] == :samplesurface
                n_samples --> Int(ceil(n_samples * scales[i]/max_scale))
            end
            o
        end
    end
end
