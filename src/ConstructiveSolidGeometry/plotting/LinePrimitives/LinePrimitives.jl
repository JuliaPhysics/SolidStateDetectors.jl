@recipe function f(ls::AbstractVector{<:AbstractLinePrimitive})
    linecolor --> :black
    xguide --> "X"
    yguide --> "Y"
    zguide --> "Z"
    aspect_ratio --> 1.0
    @series begin
        label --> "Edges"
        ls[1]
    end
    if length(ls) > 1
        for l in ls[2:end]
            @series begin
                label := nothing
                l
            end
        end
    end
end

include("Edge.jl")
include("Ellipse.jl")