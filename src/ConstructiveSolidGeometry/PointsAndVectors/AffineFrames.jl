frame_transformation(a, b, pt::AbstractCoordinatePoint) = frame_transformation(a, b)(pt)
frame_transformation(a, b, pts::AbstractVector{<:AbstractCoordinatePoint}) = frame_transformation(a, b).(pts)

frame_transformation(a, b, v::AbstractCoordinateVector) = frame_transformation(a, b)(v)
frame_transformation(a, b, vs::AbstractVector{<:AbstractCoordinateVector}) = frame_transformation(a, b).(vs)



abstract type AbstractAffineFrame end

struct GlobalAffineFrame <: AbstractAffineFrame end

const global_frame = GlobalAffineFrame()

frame_transformation(::GlobalAffineFrame, ::GlobalAffineFrame) = identity



struct LocalAffineFrame{PT<:CartesianPoint,LM} <: AbstractAffineFrame
    # origin in the global frame
    origin::PT

    # linear operator (rotation matrix, etc.) in respect to the global frame
    linop::LM
end


struct AFTransformLinOpFirst{V<:CartesianVector,LM}
    linop::LM
    offset::V
end

function (f::AFTransformLinOpFirst)(pt::CartesianPoint)
    cartesian_zero + muladd(f.linop, pt - cartesian_zero, f.offset) 
end
(f::AFTransformLinOpFirst)(v::CartesianVector) = f.linop * v

(f::AFTransformLinOpFirst)(pt::CylindricalPoint) = CylindricalPoint(f(CartesianPoint(pt)))

InverseFunctions.inverse(f::AFTransformLinOpFirst) = AFTransformOffsetFirst(-f.offset, inv(f.linop))


struct AFTransformOffsetFirst{V,M}
    offset::V
    linop::M
end

(f::AFTransformOffsetFirst)(pt::CartesianPoint) =  cartesian_zero + (f.linop * (pt - cartesian_zero + f.offset))
(f::AFTransformOffsetFirst)(v::CartesianVector) = f.linop * v

(f::AFTransformOffsetFirst)(pt::CylindricalPoint) = CylindricalPoint(f(CartesianPoint(pt)))

InverseFunctions.inverse(f::AFTransformOffsetFirst) = AFTransformLinOpFirst(inv(f.linop), -f.offset)


@inline function frame_transformation(frame::LocalAffineFrame, ::GlobalAffineFrame)
    return AFTransformLinOpFirst(frame.linop, frame.origin - cartesian_zero)
end

@inline function frame_transformation(::GlobalAffineFrame, frame::LocalAffineFrame)
    return AFTransformOffsetFirst(cartesian_zero - frame.origin, inv(frame.linop))
end

function frame_transformation(a::LocalAffineFrame, b::LocalAffineFrame)
    return frame_transformation(a, global_frame) ∘ frame_transformation(global_frame, b)
end
