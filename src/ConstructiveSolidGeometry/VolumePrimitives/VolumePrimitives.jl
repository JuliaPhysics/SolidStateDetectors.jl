ClosedPrimitive(p::VP) where {VP <:AbstractVolumePrimitive} = VP(p, COT = ClosedPrimitive, origin = p.origin, rotation = p.rotation)
OpenPrimitive(p::VP) where {VP <: AbstractVolumePrimitive} = VP(p, COT = OpenPrimitive, origin = p.origin, rotation = p.rotation)

rotation(p::AbstractVolumePrimitive) = p.rotation
origin(p::AbstractVolumePrimitive) = p.origin

distance(pt::CartesianPoint, vp::AbstractVolumePrimitive) = 
    minimum(map(p -> distance(pt, p), surfaces(vp)))

_transform_into_object_coordinate_system(pt::CartesianPoint, p::AbstractVolumePrimitive) = inv(rotation(p)) * (pt - origin(p)) 
in(pt::CartesianPoint, p::AbstractVolumePrimitive) = _in(_transform_into_object_coordinate_system(pt, p), p)
in(pt::CylindricalPoint, p::AbstractVolumePrimitive) = in(CartesianPoint(pt), p)
# Do we want to store the rotation matrix permanently in the primitive?
# We should do tests regarding the performance. It can be easily added later. 

include("Box.jl")