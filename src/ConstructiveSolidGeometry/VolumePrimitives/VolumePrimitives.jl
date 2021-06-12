include("Box.jl")

distance(pt::CartesianPoint, b::AbstractVolumePrimitive) = 
    minimum(map(p -> distance(pt, p), surfaces(b)))
