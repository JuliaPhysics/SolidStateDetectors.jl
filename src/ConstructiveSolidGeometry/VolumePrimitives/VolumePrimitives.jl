ClosedPrimitive(p::VP) where {VP <:AbstractVolumePrimitive} = VP(p, COT = ClosedPrimitive)
OpenPrimitive(p::VP) where {VP <: AbstractVolumePrimitive} = VP(p, COT = OpenPrimitive)


distance(pt::CartesianPoint, vp::AbstractVolumePrimitive) = 
minimum(map(p -> distance(pt, p), surfaces(vp)))

include("Box.jl")