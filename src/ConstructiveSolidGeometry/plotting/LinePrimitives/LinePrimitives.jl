include("Edge.jl")

@recipe function f(ls::AbstractVector{<:AbstractLinePrimitive})
    linecolor --> :black
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
