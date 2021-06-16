"""
    struct Ellipse{T,TR} <: AbstractLinePrimitive{T}

* `r::TR`: 
    * TR = Real -> Circle (a = b = r)
    * TR = (Real, Real) -> Circular Annulus (r_in = r[1], r_out = r[2])
    * TR = ((Real,), (Real,)) -> Ellipse (a = r[1][1], b = r[2][1])
    * TR = ((Real, Real),(Real, Real)) -> Elliptical Annulus \n(a_in = r[1][1], a_out = r[1][2], b_in = r[2][1], b_out = r[2][2])
    * Not all are implemented yet

* `φ::TP`: 
    * TP = Nothing <-> Full in φ
    * ...
"""
@with_kw struct Ellipse{T,TR,TP} <: AbstractLinePrimitive{T}
    r::TR = 1
    φ::TP = nothing

    origin::CartesianPoint{T} = zero(CartesianPoint{T})
    rotation::SMatrix{3,3,T,9} = one(SMatrix{3, 3, T, 9})
end

Ellipse{T,TR,TP}( e::Ellipse{T,TR,TP}; 
            origin::CartesianPoint{T} = b.origin,
            rotation::SMatrix{3,3,T,9} = b.rotation) where {T,TR,TP} =
    Ellipse{T,TR,TP}(e.r, e.φ, origin, rotation)


# function Arc(; r = 0, center = PlanarPoint(0,0), αMin = 0, αMax = 2π)
#     T = float(promote_type(typeof.((r, αMin, αMax))..., eltype(center)))
#     α = mod(T(αMax) - T(αMin), T(2π)) == 0 ? nothing : T(αMin)..T(αMax)
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
