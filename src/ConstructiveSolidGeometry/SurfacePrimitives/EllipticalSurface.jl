"""
    struct EllipticalSurface{T,TR,TP} <: AbstractPlanarSurfacePrimitive{T}

Surface primitive describing circular bases, e.g. the top or bottom base of a [`Cone`](@ref).

## Parametric types
* `T`: Precision type.
* `TR`: Type of the radius `r`.
    * `TR == T`: Full Circle (constant radius `r`, no cut-out).
    * `TR == Tuple{T, T}`: Circular Annulus (inner radius at `r[1]`, outer radius at `r[2]`).
* `TP`: Type of the angular range `φ`.
    * `TP == Nothing`: Full 2π Cone.
    * `TP == T`: Partial Elliptical Surface ranging from `0` to `φ`.
    
## Fields
* `r::TR`: Definition of the radius of the `EllipticalSurface` (in m).
* `φ::TP`: Range in polar angle `φ` over which the `EllipticalSurface` extends (in radians).
* `origin::CartesianPoint{T}`: The position of the center of the `EllipticalSurface`.
* `rotation::SMatrix{3,3,T,9}`: Matrix that describes a rotation of the `EllipticalSurface` around its `origin`.
"""
struct EllipticalSurface{T<:AbstractFloat,TR<:Union{T,Tuple{T,T}},TP<:Union{Nothing,T}} <: AbstractPlanarSurfacePrimitive{T}
    r::TR
    φ::TP 

    origin::CartesianPoint{T}
    rotation::SMatrix{3,3,T,9} 
end

#Type conversion happens here
function EllipticalSurface{T}(r, φ, origin, rotation) where {T}
    _r = _csg_convert_args(T, r)
    (_φ, _rotation) = _handle_phi(_csg_convert_args(T, φ), rotation)
    EllipticalSurface{T,typeof(_r),typeof(_φ)}(_r, _φ, origin, _rotation)
end

#Type promotion happens here
function EllipticalSurface(r::TR, φ::TP, origin::PT, rotation::ROT) where {TR<:Union{Tuple,Real}, TP, PT, ROT}
    eltypes = _csg_get_promoted_eltype.((TR, TP, PT, ROT))
    T = float(promote_type(eltypes...))
    EllipticalSurface{T}(r, φ, origin, rotation)
end

function EllipticalSurface(;
    # define default parameters as Int to not influence type promotion
    r = 1, 
    φ = nothing,
    origin = zero(CartesianPoint{Int}), 
    rotation = one(SMatrix{3, 3, Int, 9}))
    EllipticalSurface(r, φ, origin, rotation)
end

function EllipticalSurface{T}(;
    r = 1.0, 
    φ = nothing,
    origin = zero(CartesianPoint{Float64}), 
    rotation = one(SMatrix{3, 3, Float64, 9})
) where {T}
    EllipticalSurface{T}(r, φ, origin, rotation)
end

flip(es::EllipticalSurface{T,TR,Nothing}) where {T,TR} = 
    EllipticalSurface{T,TR,Nothing}(es.r, es.φ, es.origin, es.rotation * SMatrix{3,3,T,9}(1,0,0,0,-1,0,0,0,-1))
flip(es::EllipticalSurface{T,TR,T}) where {T,TR} = 
    EllipticalSurface{T,TR,T}(es.r, es.φ, es.origin, es.rotation * SMatrix{3,3,T,9}(1,0,0,0,-1,0,0,0,-1) * RotZ{T}(T(2π)-es.φ))

const CircularArea{T} = EllipticalSurface{T,T,Nothing}
const PartialCircularArea{T} = EllipticalSurface{T,T,T}

const Annulus{T} = EllipticalSurface{T,Tuple{T,T},Nothing}
const PartialAnnulus{T} = EllipticalSurface{T,Tuple{T,T},T}

Plane(es::EllipticalSurface{T}) where {T} = Plane{T}(es.origin, normal(es))
# TODO is this correct?
normal(es::EllipticalSurface{T}, ::CartesianPoint{T} = zero(CartesianPoint{T})) where {T} = 
    _transform_into_global_coordinate_system(CartesianVector{T}(zero(T), zero(T), one(T)), es)

function vertices(es::EllipticalSurface{T, T}, n_arc::Int)::Vector{CartesianPoint{T}} where {T}
    φMin, φMax = get_φ_limits(es)
    n_arc = _get_n_points_in_arc_φ(es, n_arc)
    φ = range(φMin, stop = φMax, length = n_arc+1)
    append!([es.origin],[_transform_into_global_coordinate_system(CartesianPoint{T}(es.r*cos(φ), es.r*sin(φ), 0), es) for φ in φ])
end

function vertices(es::EllipticalSurface{T, Tuple{T,T}}, n_arc::Int)::Vector{CartesianPoint{T}} where {T}
    rMin, rMax = es.r
    φMin, φMax = get_φ_limits(es)
    n_arc = _get_n_points_in_arc_φ(es, n_arc)
    φ = range(φMin, stop = φMax, length = n_arc+1)
    [_transform_into_global_coordinate_system(CartesianPoint{T}(r*cos(φ), r*sin(φ), 0), es) for r in (rMin,rMax) for φ in φ]
end

function sample(es::EllipticalSurface{T}, spacing::T)::Vector{CartesianPoint{T}} where {T}
    φMin, φMax = get_φ_limits(es)
    Δφ = abs(φMax - φMin)
    full2π = isnothing(es.φ)
    rMin, rMax = length(es.r) == 1 ? (0, es.r) : es.r
    [_transform_into_global_coordinate_system(CartesianPoint{T}(r*cos(φ), r*sin(φ), 0), es) 
        for r in range(rMin, stop = rMax, length = max(2,1+Int(ceil((rMax-rMin)/spacing)))) 
        for φ in (r == 0 ? [φMin] : range(φMin, stop = φMax - full2π*spacing/r, length = max(2,1+Int(ceil(Δφ*r/spacing)))))]
end

function connections(es::EllipticalSurface{T, T}, n_arc::Int)::Vector{Vector{Int}} where {T}
    n_arc = _get_n_points_in_arc_φ(es, n_arc)
    [[1,i,i+1] for i in 2:n_arc+1]
end

function connections(es::EllipticalSurface{T, Tuple{T,T}}, n_arc::Int)::Vector{Vector{Int}} where {T}
    n_arc = _get_n_points_in_arc_φ(es, n_arc)
    [[i,i+1,i+n_arc+2,i+n_arc+1] for i in 1:n_arc]
end

function connections(es::EllipticalSurface{T, T}, n_arc::Int, n_vert_lines::Int)::Vector{Vector{Int}} where {T}
    n_arc = _get_n_points_in_arc_φ(es, n_arc)
    radii = [[1, i + 1] for i in _get_vert_lines_range(es,n_arc,n_vert_lines)]
    circ =  [[i, i+1] for i in 2:n_arc+1]
    append!(radii, circ)
end

function connections(es::EllipticalSurface{T, Tuple{T,T}}, n_arc::Int, n_vert_lines::Int)::Vector{Vector{Int}} where {T}
    n_arc = _get_n_points_in_arc_φ(es, n_arc)
    radii = [[i, i + n_arc + 1] for i in _get_vert_lines_range(es,n_arc,n_vert_lines)]
    circ1 =  [[i, i + 1] for i in 1:n_arc]
    circ2 =  [[i, i + 1] for i in n_arc+2:2*n_arc+1]
    append!(radii, circ1, circ2)
end

get_label_name(::EllipticalSurface) = "Elliptical Surface"

extremum(es::EllipticalSurface{T,T}) where {T} = es.r
extremum(es::EllipticalSurface{T,Tuple{T,T}}) where {T} = es.r[2] # r_out always larger r_in: es.r[2] > es.r[2]

get_φ_limits(es::EllipticalSurface{T,<:Any,T}) where {T} = T(0), es.φ
get_φ_limits(cm::EllipticalSurface{T,<:Any,Nothing}) where {T} = T(0), T(2π)

function lines(sp::CircularArea{T}; n = 2) where {T} 
    circ = Circle{T}(sp.r, sp.φ, sp.origin, sp.rotation)
    φs = range(T(0), step = T(2π) / n, length = n)
    edges = [ Edge{T}(_transform_into_global_coordinate_system(CartesianPoint(CylindricalPoint{T}(zero(T), φ, zero(T))), sp),
                      _transform_into_global_coordinate_system(CartesianPoint(CylindricalPoint{T}(sp.r, φ, zero(T))), sp)) for φ in φs ]
    return (circ, edges)
end
function lines(sp::PartialCircularArea{T}; n = 2) where {T} 
    circ = PartialCircle{T}(sp.r, sp.φ, sp.origin, sp.rotation)
    φs = range(T(0), stop = sp.φ, length = n)
    edges = [ Edge{T}(_transform_into_global_coordinate_system(CartesianPoint(CylindricalPoint{T}(zero(T), φ, zero(T))), sp),
                      _transform_into_global_coordinate_system(CartesianPoint(CylindricalPoint{T}(sp.r, φ, zero(T))), sp)) for φ in φs ]
    return (circ, edges)
end

function lines(sp::Annulus{T}; n = 2) where {T} 
    circ_in  = Circle{T}(sp.r[1], sp.φ, sp.origin, sp.rotation)
    circ_out = Circle{T}(sp.r[2], sp.φ, sp.origin, sp.rotation)
    φs = range(T(0), step = T(2π) / n, length = n)
    edges = [ Edge{T}(_transform_into_global_coordinate_system(CartesianPoint(CylindricalPoint{T}(sp.r[1], φ, zero(T))), sp),
                      _transform_into_global_coordinate_system(CartesianPoint(CylindricalPoint{T}(sp.r[2], φ, zero(T))), sp)) for φ in φs ]
    return (circ_in, circ_out, edges)
end
function lines(sp::PartialAnnulus{T}; n = 2) where {T} 
    circ_in  = PartialCircle{T}(sp.r[1], sp.φ, sp.origin, sp.rotation)
    circ_out = PartialCircle{T}(sp.r[2], sp.φ, sp.origin, sp.rotation)
    φs = range(T(0), stop = sp.φ, length = n)
    edges = [ Edge{T}(_transform_into_global_coordinate_system(CartesianPoint(CylindricalPoint{T}(sp.r[1], φ, zero(T))), sp),
                      _transform_into_global_coordinate_system(CartesianPoint(CylindricalPoint{T}(sp.r[2], φ, zero(T))), sp)) for φ in φs ]
    return (circ_in, circ_out, edges)
end

function distance_to_surface(pt::AbstractCoordinatePoint{T}, a::EllipticalSurface{T, <:Any, TP})::T where {T, TP<:Union{Nothing,T}}

    pt_cart = _transform_into_object_coordinate_system(CartesianPoint(pt), a)
    pt_cyl = CylindricalPoint(pt_cart)
    rMin, rMax = _radial_endpoints(a.r)

    # Full ellipse
    return if isnothing(a.φ) || _in_φ(pt_cyl, a.φ)
        hypot(pt_cyl.z, max(zero(T), rMin - pt_cyl.r, pt_cyl.r - rMax))
    else
        # Nearest φ-boundary is either 0 or a.φ
        φNear = _φNear(pt_cyl.φ, a.φ)
        p1 = CartesianPoint{T}(rMin*cos(φNear), rMin*sin(φNear), 0)
        p2 = CartesianPoint{T}(rMax*cos(φNear), rMax*sin(φNear), 0)
        distance_to_line(pt_cart, Edge(p1, p2))
    end
end
