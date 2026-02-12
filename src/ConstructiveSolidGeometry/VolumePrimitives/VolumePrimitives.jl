# All volume primitives need to have the following methods defined: 
# * _in(::CartesianPoint, ::VolumePrimitive)  (in the local coordinate system)
# * surfaces(::VolumePrimitive) --> all surfaces of the volume primitive as Tuple
# * Geometry(...) --> a method to read the geometry from a config file
# * Dictionary(...) --> a method to convert the primitive to a config 

ClosedPrimitive(p::VP) where {VP <:AbstractVolumePrimitive} = VP(p, COT = ClosedPrimitive)
OpenPrimitive(p::VP) where {VP <: AbstractVolumePrimitive} = VP(p, COT = OpenPrimitive)
switchClosedOpen(p::AbstractVolumePrimitive{<:Any,ClosedPrimitive}) = OpenPrimitive(p)
switchClosedOpen(p::AbstractVolumePrimitive{<:Any,OpenPrimitive}) = ClosedPrimitive(p)

isClosedPrimitive(vp::AbstractVolumePrimitive{<:Any,ClosedPrimitive}) = true
isClosedPrimitive(vp::AbstractVolumePrimitive{<:Any,OpenPrimitive}) = false

distance(pt::CartesianPoint, vp::AbstractVolumePrimitive) = 
minimum(map(p -> distance(pt, p), surfaces(vp)))

function sample(vp::AbstractVolumePrimitive) 
    surfs = surfaces(vp)
    ls = unique!(collect(Base.Iterators.flatten(lines.(surfs))))
    unique!(mapreduce(sample, vcat, ls))
end

function surfaces(vp::AbstractVolumePrimitive{<:Any,OpenPrimitive})
    flip.(surfaces(ClosedPrimitive(vp)))
end

include("Box.jl")
include("Cone.jl")
include("Ellipsoid.jl")
include("RegularPrism.jl")
include("Polycone.jl")
include("Torus.jl")