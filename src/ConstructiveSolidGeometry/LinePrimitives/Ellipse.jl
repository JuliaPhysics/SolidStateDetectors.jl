# """
#     struct Ellipse{T,TR} <: AbstractLinePrimitive{T}
# 
# * `r::TR`: 
#     * TR = Real -> Circle (a = b = r)
#     * TR = (Real, Real) -> Circular Annulus (r_in = r[1], r_out = r[2])
#     * TR = ((Real,), (Real,)) -> Ellipse (a = r[1][1], b = r[2][1])
#     * TR = ((Real, Real),(Real, Real)) -> Elliptical Annulus \n(a_in = r[1][1], a_out = r[1][2], b_in = r[2][1], b_out = r[2][2])
#     * Not all are implemented yet
# 
# * `φ::TP`: 
#     * TP = Nothing <-> Full in φ
#     * ...
# """
struct Ellipse{T,TR,TP<:Union{Nothing,T}} <: AbstractLinePrimitive{T}
    r::TR 
    φ::TP 

    origin::CartesianPoint{T}
    rotation::SMatrix{3,3,T,9}
end

#Type conversion happens here
function Ellipse{T}(r::TR, φ, origin, rotation) where {T,TR<:Union{Real,Tuple}}
    _r = _csg_convert_args(T, r)
    _φ = _csg_convert_args(T, φ)
    Ellipse{T,typeof(_r),typeof(_φ)}(_r, _φ, origin, rotation)
end

#Type promotion happens here
function Ellipse(r::TR, φ::TP, origin::PT, rotation::ROT) where {TR<:Union{Real,Tuple}, TP<:Union{Nothing,Real}, PT<:(CartesianPoint), ROT<:(StaticArrays.SMatrix{3, 3, T, 9} where T)}
    eltypes = _csg_get_promoted_eltype.((TR, TP, PT, ROT))
    T = float(promote_type(eltypes...))
    Ellipse{T}(r, φ, origin, rotation)
end

function Ellipse(;
    r = 1,
    φ = nothing,
    origin = zero(CartesianPoint{Int}), 
    rotation = one(SMatrix{3, 3, Int, 9})
) 
    Ellipse(r, φ, origin, rotation)
end

function Ellipse{T}(;
    r = 1.0,
    φ = nothing,
    origin = zero(CartesianPoint{Float64}), 
    rotation = one(SMatrix{3, 3, Float64, 9})
) where {T}
    Ellipse{T}(r, φ, origin, rotation)
end

const Circle{T} = Ellipse{T,T,Nothing}
const PartialCircle{T} = Ellipse{T,T,T}

extremum(e::Ellipse{T,T}) where {T} = e.r
extremum(e::Ellipse{T,Tuple{T,T}}) where {T} = max(e.r...)

function sample(e::Circle{T}; n = 4)::Vector{CartesianPoint{T}} where {T}
    φs = range(T(0), step = T(2π) / n, length = n)
    pts = Vector{CartesianPoint{T}}(undef, n)
    for i in eachindex(pts)
        pts[i] = _transform_into_global_coordinate_system(CartesianPoint(CylindricalPoint{T}(e.r, φs[i], zero(T))), e)
    end
    pts
end
function sample(e::PartialCircle{T}; n = 2)::Vector{CartesianPoint{T}} where {T}
    φs = range(T(0), stop = e.φ, length = n)
    pts = Vector{CartesianPoint{T}}(undef, n)
    for i in eachindex(pts)
        pts[i] = _transform_into_global_coordinate_system(CartesianPoint(CylindricalPoint{T}(e.r, φs[i], zero(T))), e)
    end
    pts
end



# function Arc(; r = 0, center = PlanarPoint(0,0), αMin = 0, αMax = 2π)
#     T = float(promote_type(typeof.((r, αMin, αMax))..., eltype(center)))
#     α = mod(T(αMax    ) - T(αMin), T(2π)) == 0 ? nothing : T(αMin)..T(αMax)
#     Arc(T, T(r), PlanarPoint{T}(center), α)
# end

# Arc(r, center, αMin, αMax) = Arc(; r = r, center = center, αMin = αMin, αMax = αMax)

# Circle(r::T, center::PlanarPoint{T}) where {T} = Arc(T, r, center, nothing)
# Circle(a::Arc{T}) where {T} = Arc(T, a.r, a.center, nothing)

# get_α_at_u_v(a::Arc{T}, u::Real, v::Real) where {T} = mod(atan(v - a.center.v, u - a.center.u), 2π) #u,v are planar coordinates

# get_α_limits(a::Arc{T, Nothing}) where {T} = (T(0), T(2π), true)
# get_α_limits(a::Arc{T, <:AbstractInterval}) where {T} = (a.α.left, a.α.right, false)

# distance_to_line(point::PlanarPoint{T}, a::Arc{T, Nothing}) where {T} = abs(norm(point - a.center) - a.r)

# function distance_to_line(point::PlanarPoint{T}, a::Arc{T, <:AbstractInterval})::T where {T}
#     αMin::T, αMax::T, _ = get_α_limits(a)
#     p_t = point - a.center
#     α = atan(p_t.v, p_t.u)
#     if _in_angular_interval_closed(α, a.α)
#         return abs(norm(p_t) - a.r)
#     else
#         sαNear, cαNear = sincos(_φNear(α, αMin, αMax))
#         return norm(p_t - PlanarPoint{T}(a.r*cαNear, a.r*sαNear))
#     end
# end
