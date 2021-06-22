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

@recipe function f(contact::Contact)
    linecolor --> contact.id
    l = contact.name != "" ? "$(contact.name) (id: $(contact.id))" : "Contact - id: $(contact.id)"
    label --> l
    contact.geometry
end

@recipe function f(det::SolidStateDetector)
    xguide --> "x / m"
    yguide --> "y / m"
    zguide --> "z / m"
    show_normal --> false

    @series begin
        det.semiconductor
    end
    for c in det.contacts
        @series begin
            c
        end
    end
    if !ismissing(det.passives)
        for p in det.passives
            @series begin
                p
            end
        end
    end
end
