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


struct AffineFrameTransformation{V<:CartesianVector,LM}
    offset::V
    linop::LM
end

(f::AffineFrameTransformation)(pt::CartesianPoint) = muladd(f.linop, pt, f.offset)
(f::AffineFrameTransformation)(v::CartesianVector) = f.linop * v

(f::AffineFrameTransformation)(pt::CylindricalPoint) = CylindricalPoint(f(CartesianPoint(pt)))

InverseFunctions.inverse(f::AffineFrameTransformation) = InvAffineFrameTransformation(f.offset, f.linop)


struct InvAffineFrameTransformation{V,M}
    offset::V
    linop::M
end

(f::InvAffineFrameTransformation)(pt::CartesianPoint) = f.linop \ (pt - f.offset)
(f::InvAffineFrameTransformation)(v::CartesianVector) = f.linop \ v

(f::InvAffineFrameTransformation)(pt::CylindricalPoint) = CylindricalPoint(f(CartesianPoint(pt)))

InverseFunctions.inverse(f::InvAffineFrameTransformation) = AffineFrameTransformation(f.offset, f.linop)


@inline function frame_transformation(frame::LocalAffineFrame, ::GlobalAffineFrame)
    return AffineFrameTransformation(frame.origin - cartesian_zero, frame.linop)
end

@inline function frame_transformation(::GlobalAffineFrame, frame::LocalAffineFrame)
    return InvAffineFrameTransformation(frame.origin - cartesian_zero, frame.linop)
end

function frame_transformation(a::LocalAffineFrame, b::LocalAffineFrame)
    return frame_transformation(a, global_frame) ∘ frame_transformation(global_frame, b)
end
