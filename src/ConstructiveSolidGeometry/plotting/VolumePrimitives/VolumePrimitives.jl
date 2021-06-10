@recipe function f(p::AbstractVolumePrimitive)
    fs = faces(p)
    # In principle it would be better to directly get the edges of the primitive
    # Can be improved later...
    linecolor --> :black
    @series begin
        label --> "$(nameof(typeof(p)))"
        fs
    end
end

@recipe function f(p1::AbstractVolumePrimitive, p2::AbstractVolumePrimitive)
    fs1 = faces(p1)
    fs2 = faces(p2)
    # In principle it would be better to directly get the edges of the primitive
    # Can be improved later...
    linecolor --> :black
    @series begin
        if plotattributes[:csguniontype] == :intersection
            linestyle := :dot
        else
            linestyle := :solid
        end
        fs1
    end
    @series begin
        if plotattributes[:csguniontype] == :intersection
            linestyle := :dot
        elseif plotattributes[:csguniontype] == :difference
            linestyle := :dash
        else
            linestyle := :solid
        end
        fs2
    end
end