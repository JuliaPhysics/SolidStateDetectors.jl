"""
    struct EllipsoidMantle{T,TR,TP,TT,D} <: AbstractCurvedSurfacePrimitive{T}

Surface primitive describing the surface of an [`Ellipsoid`](@ref).

## Parametric types
* `T`: Precision type.
* `TR`: Type of the radius `r`.
    * `TR == T`: SphereMantle (constant radius `r` along all axes).
* `TP`: Type of the azimuthial angle `φ`.
    * `TP == Nothing`: Full 2π in `φ`.
* `TT`: Type of the polar angle `θ`.
    * `TT == Nothing`: Full 2π in `θ`.
* `D`: Direction in which the normal vector points (`:inwards` or `:outwards`).

## Fields
* `r::TR`: Definition of the radius of the `EllipsoidMantle` (in m).
* `φ::TP`: Range in azimuthial angle `φ` of the `EllipsoidMantle`.
* `θ::TT`: Range in polar angle `θ` of the `EllipsoidMantle`.
* `origin::CartesianPoint{T}`: Origin of the `Ellipsoid` which has this `EllipsoidMantle` as surface.
* `rotation::SMatrix{3,3,T,9}`: Rotation matrix of the `Ellipsoid` which has this `EllipsoidMantle` as surface.
"""
struct EllipsoidMantle{T,TR,TP,TT,D} <: AbstractCurvedSurfacePrimitive{T}
    r::TR
    φ::TP
    θ::TT

    origin::CartesianPoint{T}
    rotation::SMatrix{3,3,T,9}
end

#Type conversion happens here
function EllipsoidMantle{T}(D, r, φ::Nothing, θ::Nothing, origin, rotation) where {T}
    _r = _csg_convert_args(T, r)
    EllipsoidMantle{T,typeof(_r),Nothing,Nothing,D}(_r, φ, θ, origin, rotation)
end

#Type promotion happens here
function EllipsoidMantle(D,r::TR, φ::TP, θ::TT, origin::PT, rotation::ROT) where {TR, TP, TT, PT, ROT}
    eltypes = _csg_get_promoted_eltype.((TR, TP, TT, PT, ROT))
    T = float(promote_type(eltypes...))
    EllipsoidMantle{T}(D, r, φ, θ, origin, rotation)
end

function EllipsoidMantle(D::Symbol=:outwards;
    # define default parameters as Int to not influence type promotion
    r = 1, 
    φ = nothing,
    θ = nothing,
    origin = zero(CartesianPoint{Int}), 
    rotation = one(SMatrix{3, 3, Int, 9}))
    EllipsoidMantle(D, r, φ, θ, origin, rotation)
end

function EllipsoidMantle{T}(D::Symbol=:outwards;
    r = 1.0, 
    φ = nothing,
    θ = nothing,
    origin = zero(CartesianPoint{Float64}), 
    rotation = one(SMatrix{3, 3, Float64, 9})
) where {T}
    EllipsoidMantle{T}(D, r, φ, θ, origin, rotation)
end

flip(em::EllipsoidMantle{T,TR,TP,TT,:inwards}) where {T,TR,TP,TT} = 
    EllipsoidMantle{T,TR,TP,TT,:outwards}(em.r, em.φ, em.θ, em.origin, em.rotation )

const FullSphereMantle{T,D} = EllipsoidMantle{T,T,Nothing,Nothing,D}
const FullEllipsoidMantle{T,D} = EllipsoidMantle{T,NTuple{3,T},Nothing,Nothing,D}

get_φ_limits(em::EllipsoidMantle{T,<:Any,Tuple{T,T}}) where {T} = em.φ[1], em.φ[2]
get_φ_limits(em::EllipsoidMantle{T,<:Any,Nothing}) where {T} = T(0), T(2π)

get_θ_limits(em::EllipsoidMantle{T,<:Any,<:Any,Tuple{T,T}}) where {T} = em.θ[1], em.θ[2]
get_θ_limits(em::EllipsoidMantle{T,<:Any,<:Any,Nothing}) where {T} = T(-π/2), T(π/2)

get_radii(em::EllipsoidMantle{T,T}) where {T} = (em.r, em.r, em.r)
get_radii(em::EllipsoidMantle{T,NTuple{3,T}}) where {T} = em.r

function lines(em::FullSphereMantle{T}) where {T} 
    ellipse_xy = Ellipse{T,T,Nothing}(em.r[1], em.φ, em.origin, em.rotation)
    ellipse_xz = Ellipse{T,T,Nothing}(em.r[1], em.φ, em.origin, em.rotation * RotX(T(π)/2))
    ellipse_yz = Ellipse{T,T,Nothing}(em.r[1], em.φ, em.origin, em.rotation * RotX(T(π)/2) * RotY(T(π)/2))
    (ellipse_xy, ellipse_xz, ellipse_yz)
end
function lines(em::FullEllipsoidMantle{T}) where {T} 
    ellipse_xy = Ellipse{T,NTuple{2, Tuple{T}},Nothing}(((em.r[1],), (em.r[2],)), em.φ, em.origin, em.rotation)
    ellipse_xz = Ellipse{T,NTuple{2, Tuple{T}},Nothing}(((em.r[1],), (em.r[3],)), em.φ, em.origin, em.rotation * RotX(T(π)/2))
    ellipse_yz = Ellipse{T,NTuple{2, Tuple{T}},Nothing}(((em.r[2],), (em.r[3],)), em.φ, em.origin, em.rotation * RotX(T(π)/2) * RotY(T(π)/2))
    (ellipse_xy, ellipse_xz, ellipse_yz)
end

extremum(e::EllipsoidMantle{T,T}) where {T} = e.r
extremum(e::EllipsoidMantle{T,NTuple{3,T}}) where {T} = max(e.r...)

function normal(em::EllipsoidMantle{T,NTuple{3,T},TP,TT,:outwards}, pt::CartesianPoint{T})::CartesianVector{T} where {T,TP,TT}
    # not normalized, do we want this?
    # Or wrap this into somehting like `normal(em, pt) = normalize(direction(em, pt))` ?
    p = _transform_into_object_coordinate_system(pt, em)
    obj_normal = CartesianVector{T}(sign(p.x)*(p.x/em.r[1])^2, sign(p.y)*(p.y/em.r[2])^2, sign(p.z)*(p.z/em.r[3])^2) # We might want to store the inv(em.r) in the struct?
    _transform_into_global_coordinate_system(obj_normal, em)
end
function normal(em::EllipsoidMantle{T,T,TP,TT,:outwards}, pt::CartesianPoint{T})::CartesianVector{T} where {T,TP,TT}
    # not normalized, do we want this?
    # Or wrap this into somehting like `normal(em, pt) = normalize(direction(em, pt))` ?
    p = _transform_into_object_coordinate_system(pt, em)
    obj_normal = CartesianVector{T}(sign(p.x)*(p.x/em.r)^2, sign(p.y)*(p.y/em.r)^2, sign(p.z)*(p.z/em.r)^2) # We might want to store the inv(em.r) in the struct?
    _transform_into_global_coordinate_system(obj_normal, em)
end
normal(em::EllipsoidMantle{T,TR,TP,TT,:inwards}, pt::CartesianPoint{T}) where {T,TR,TP,TT} = -normal(flip(em), pt)

function vertices(em::EllipsoidMantle{T}, n_arc::Int)::Vector{CartesianPoint{T}} where {T}
    rx, ry, rz = get_radii(em) 
    φMin, φMax = get_φ_limits(em)
    θMin, θMax = get_θ_limits(em)
    n_arcφ = _get_n_points_in_arc_φ(em, n_arc) 
    n_arcθ = _get_n_points_in_arc_θ(em, n_arc)
    
    θ = range(θMin, stop = θMax, length = n_arcθ + 1)
    φ = range(φMin, stop = φMax, length = n_arcφ + 1)
    
    [_transform_into_global_coordinate_system(CartesianPoint{T}(rx*cos(θ)*cos(φ), ry*cos(θ)*sin(φ), rz*sin(θ)), em) for θ in θ for φ in φ]
end

function sample(em::EllipsoidMantle{T}, spacing::T)::Vector{CartesianPoint{T}} where {T}
    rx, ry, rz = get_radii(em) 
    r = (rx*ry*rz)^(1/3)
    φMin, φMax = get_φ_limits(em)
    θMin, θMax = get_θ_limits(em)
    Δφ = abs(φMax - φMin)
    Δθ = abs(θMax - θMin)
    full2π = isnothing(em.φ)
    
    [_transform_into_global_coordinate_system(CartesianPoint{T}(rx*cos(θ)*cos(φ), ry*cos(θ)*sin(φ), rz*sin(θ)), em) 
        for θ in range(θMin, stop = θMax, length = max(2,1+Int(ceil(Δθ*r/spacing)))) 
        for φ in (mod(θ, T(2π)) in T[π/2,3π/2] ? [φMin] : range(φMin, stop = φMax - full2π*spacing/(r*cos(θ)), length = max(2,1+Int(ceil(Δφ*r*cos(θ)/spacing)))))]
end

function connections(em::EllipsoidMantle, n_arc::Int)::Vector{Vector{Int}}
    n_arcφ = _get_n_points_in_arc_φ(em, n_arc) 
    n_arcθ = _get_n_points_in_arc_θ(em, n_arc)
    [[i+(n_arcφ+1)*j,i+1+(n_arcφ+1)*j,i+1+(n_arcφ+1)*(j+1),i+(n_arcφ+1)*(j+1)] for j in 0:n_arcθ-1 for i in 1:n_arcφ]
end

function _get_hor_lines_idx(em::EllipsoidMantle{T}, n_arcθ::Int)::Vector{Int} where {T}
    θMin, θMax = get_θ_limits(em) 
    θMin, θMax = minmax(mod(θMin, T(2π)), mod(θMax, T(2π)))
    if θMin == T(π/2) && θMax == T(3π/2) 
        [Int(floor(n_arcθ/2)) + 1]
    else
        [1, n_arcθ+1]
    end
end

function connections(em::EllipsoidMantle, n_arc::Int, n_vert_lines::Int)::Vector{Vector{Int}} 
    n_arcφ = _get_n_points_in_arc_φ(em, n_arc) 
    n_arcθ = _get_n_points_in_arc_θ(em, n_arc)
    vcircs = [[i*(n_arcφ+1)+c, (i+1)*(n_arcφ+1)+c] for c in _get_vert_lines_range(em,n_arcφ,n_vert_lines) for i in 0:n_arcθ-1]
    hcircs =  [[i + (c-1)*(n_arcφ+1), i + (c-1)*(n_arcφ+1) + 1] for c in _get_hor_lines_idx(em, n_arcθ) for i in 1:n_arcφ]
    append!(vcircs, hcircs)
end

get_label_name(::EllipsoidMantle) = "Ellipsoid Mantle"

"""
    intersection(em::EllipsoidMantle{T}, l::Line{T}) where {T}

Calculates the intersections of a `Line` with a `EllipsoidMantle`.

## Arguments
* `cm::EllipsoidMantle{T}`: The `EllipsoidMantle`.
* `l::Line{T}`: The `Line`.

!!! note 
    The function will always return 2 CartesianPoint's.
    If the line just touches the mantle, the two points will be the same. 
    If the line does not touch the mantle at all, the two points will have NaN's as there coordinates.
"""
function intersection(em::EllipsoidMantle{T,NTuple{3,T}}, l::Line{T}) where {T}
    obj_l = _transform_into_object_coordinate_system(l, em) # direction is not normalized
    
    L1 = obj_l.origin.x
    L2 = obj_l.origin.y
    L3 = obj_l.origin.z
    D1 = obj_l.direction.x
    D2 = obj_l.direction.y
    D3 = obj_l.direction.z
    
    R1 = em.r[1]
    R2 = em.r[2]
    R3 = em.r[3]

    term3 = D1^2/R1^2 + D2^2/R2^2 + D3^2/R3^2
    term1 = ((2D1*L1)/R1^2 + (2D2*L2)/R2^2 + (2D3*L3)/R3^2)^2 - 4*term3*(L1^2/R1^2 + L2^2/R2^2 + L3^2/R3^2 - 1)
    term2 = -(2D1*L1)/R1^2 - (2D2*L2)/R2^2 - (2D3*L3)/R3^2

    term1 = term1 < 0 ? sqrt(abs(term1)) : sqrt(term1)

    λ1 = (term2 - term1) / 2term3
    λ2 = (term2 + term1) / 2term3

    ints1 = obj_l.origin + λ1 * obj_l.direction 
    ints2 = obj_l.origin + λ2 * obj_l.direction 
    return _transform_into_global_coordinate_system(ints1, em), 
           _transform_into_global_coordinate_system(ints2, em)
end

function intersection(em::EllipsoidMantle{T,T}, l::Line{T}) where {T}
    obj_l = _transform_into_object_coordinate_system(l, em) # direction is not normalized
    
    L1 = obj_l.origin.x
    L2 = obj_l.origin.y
    L3 = obj_l.origin.z
    D1 = obj_l.direction.x
    D2 = obj_l.direction.y
    D3 = obj_l.direction.z
    
    R = em.r

    term3 = D1^2 + D2^2 + D3^2
    term1 = (2D1*L1 + 2D2*L2 + 2D3*L3)^2 - 4*term3*(L1^2 + L2^2 + L3^2 - R^2)
    term2 = -D1*L1 - D2*L2 - D3*L3
    
    if term1 < 0 term1 = abs(term1) end

    λ1 = (term2 - sqrt(term1)/2) / term3
    λ2 = (term2 + sqrt(term1)/2) / term3

    ints1 = obj_l.origin + λ1 * obj_l.direction 
    ints2 = obj_l.origin + λ2 * obj_l.direction 
    return _transform_into_global_coordinate_system(ints1, em), 
           _transform_into_global_coordinate_system(ints2, em)
end


function distance_to_surface(pt::AbstractCoordinatePoint{T}, s::FullSphereMantle{T})::T where {T}
    abs(norm(_transform_into_object_coordinate_system(CartesianPoint(pt), s) - cartesian_zero) - s.r)
end
