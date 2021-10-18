@recipe function f(p::Passive)
    linecolor --> :grey
    l = p.name != "" ? p.name : "Passive $(p.id)"
    label --> l
    p.geometry
end
@recipe function f(sc::Semiconductor)
    linecolor --> :black
    label --> "Semiconductor"
    sc.geometry
end

@recipe function f(contact::Contact{T}) where {T}
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

@recipe function f(det::SolidStateDetector; show_semiconductor = false, show_passives = true)

    show_normal --> false

    if show_semiconductor
        @series begin
            det.semiconductor
        end
    end
    for c in det.contacts
        @series begin
            c
        end
    end
    if show_passives && !ismissing(det.passives)
        for p in det.passives
            @series begin
                p
            end
        end
    end
end
