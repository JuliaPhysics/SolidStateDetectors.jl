@recipe function f(p::AbstractVolumePrimitive)
    fs = surfaces(p)
    # In principle it would be better to directly get the edges of the primitive
    # Can be improved later...
    linecolor --> :black
    @series begin
        label --> "$(nameof(typeof(p)))"
        [fs...]
    end
end

