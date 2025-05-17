abstract type AbstractAffineFrame end

affine_frame(frame::AbstractAffineFrame) = frame

struct GlobalAffineFrame <: AbstractAffineFrame end

const global_frame = GlobalAffineFrame()

frame_transformation(a, b, pt::CartesianPoint) = frame_transformation(a, b)(pt)
frame_transformation(a, b, pt::AbstractVector{<:CartesianPoint}) = frame_transformation(a, b).(pt)

frame_transformation(::GlobalAffineFrame, ::GlobalAffineFrame) = identity



struct LocalAffineFrame{T} <: AbstractAffineFrame
    global_origin::CartesianPoint{T}
    global_rotation::CartesianRotation{T}
end


struct AffineFrameTransformation{V,M}
    translation::V
    rotation::M
end

(f::AffineFrameTransformation)(pt::CartesianPoint) = muladd(f.rotation, pt, f.translation)
(f::AffineFrameTransformation)(v::CartesianVector) = f.rotation * v

InverseFunctions.inverse(f::AffineFrameTransformation) = InvAffineFrameTransformation(f.translation, f.rotation)

@inline function frame_transformation(frame::LocalAffineFrame, ::GlobalAffineFrame)
    pt = frame.global_origin
    # semantically `v = pt - CartesianVector(0, 0, 0)`, but faster:
    trans = CartesianVector(pt.x, pt.y, pt.z)
    rot = frame.global_rotation
    return AffineFrameTransformation(trans, rot)
end



struct InvAffineFrameTransformation{V,M}
    translation::V
    rotation::M
end

(f::InvAffineFrameTransformation)(pt::CartesianPoint) = f.rotation \ (pt - f.translation)
(f::InvAffineFrameTransformation)(v::CartesianVector) = f.rotation \ v

InverseFunctions.inverse(f::InvAffineFrameTransformation) = AffineFrameTransformation(f.translation, f.rotation)

@inline function frame_transformation(::GlobalAffineFrame, frame::LocalAffineFrame)
    return inverse(frame_transformation(global_frame, frame))
end
