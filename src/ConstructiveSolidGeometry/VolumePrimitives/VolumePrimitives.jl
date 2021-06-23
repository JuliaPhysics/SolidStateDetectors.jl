ClosedPrimitive(p::VP) where {VP <:AbstractVolumePrimitive} = VP(p, COT = ClosedPrimitive)
OpenPrimitive(p::VP) where {VP <: AbstractVolumePrimitive} = VP(p, COT = OpenPrimitive)
switchClosedOpen(p::AbstractVolumePrimitive{<:Any,ClosedPrimitive}) = OpenPrimitive(p)
switchClosedOpen(p::AbstractVolumePrimitive{<:Any,OpenPrimitive}) = ClosedPrimitive(p)

distance(pt::CartesianPoint, vp::AbstractVolumePrimitive) = 
minimum(map(p -> distance(pt, p), surfaces(vp)))

function sample(vp::AbstractVolumePrimitive) 
    surfs = surfaces(vp)
    ls = unique!(collect(Base.Iterators.flatten(lines.(surfs))))
    unique!(mapreduce(sample, vcat, ls))
end

include("Box.jl")
include("Cone.jl")
include("RegularPrism.jl")