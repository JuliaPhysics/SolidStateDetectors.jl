@recipe function f(p::AbstractVolumePrimitive)
    fs = surfaces(p)
    # In principle it would be better to directly get the edges of the primitive
    # Can be improved later...
    linecolor --> :black
    xguide --> "X"
    yguide --> "Y"
    zguide --> "Z"
    aspect_ratio --> 1.0
    @series begin
        label --> "$(nameof(typeof(p)))"
        [fs...]
    end
end

