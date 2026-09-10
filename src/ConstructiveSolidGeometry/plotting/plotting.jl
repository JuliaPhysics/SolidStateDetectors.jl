get_label_name(p::AbstractPrimitive) = string(nameof(typeof(p)))

include("PointsAndVectors/PointsAndVectors.jl")
include("LinePrimitives/LinePrimitives.jl")
include("SurfacePrimitives/SurfacePrimitives.jl")
include("VolumePrimitives/VolumePrimitives.jl")
include("Meshing.jl")
include("CSG/CSG.jl")
