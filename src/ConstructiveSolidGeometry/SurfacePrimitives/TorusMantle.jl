"""
    struct TorusMantle{T,TP,TT,D} <: AbstractCurvedSurfacePrimitive{T}

Surface primitive describing the mantle of a [`Torus`](@ref).

## Parametric types
* `T`: Precision type.
* `TP`: Type of the azimuthial angle `φ`.
    * `TP == Nothing`: Full 2π in `φ`.
    * `TP == T`: Partial Torus Mantle ranging from `0` to `φ`.
* `TT`: Type of the polar angle `θ`.
    * `TT == Nothing`: Full 2π in `θ`.
    * `TP == Tuple{T,T}`: Partial Torus Mantle ranging from `θ[1]` to `θ[2]`.
* `D`: Direction in which the normal vector points (`:inwards` or `:outwards`).

    
## Fields
* `r_torus::T`: Distance of the center of the `TorusMantle` to the center of the tube (in m).
* `r_tube::T`: Radius of the tube of the `TorusMantle` (in m).
* `φ::TP`: Range in azimuthial angle `φ` of the `TorusMantle`.
* `θ::TT`: Range in polar angle `θ` of the `TorusMantle`.
* `origin::CartesianPoint{T}`: The position of the center of the `TorusMantle`.
* `rotation::SMatrix{3,3,T,9}`: Matrix that describes a rotation of the `TorusMantle` around its `origin`.
"""
struct TorusMantle{T,TP<:Union{Nothing,T},TT,D} <: AbstractCurvedSurfacePrimitive{T}
    r_torus::T 
    r_tube::T 
    φ::TP 
    θ::TT 

    origin::CartesianPoint{T} 
    rotation::SMatrix{3,3,T,9} 
end

#Type conversion happens here
function TorusMantle{T}(D,r_torus, r_tube, φ::TP, θ, origin, rotation) where {T,TP}
    _r_torus = _csg_convert_args(T, r_torus)
    _r_tube = _csg_convert_args(T, r_tube)
    (_φ, _rotation) = _handle_phi(_csg_convert_args(T, φ), rotation)
    _θ = _csg_convert_args(T, θ)
    TorusMantle{T,typeof(_φ),typeof(_θ),D}(_r_torus, _r_tube, _φ, _θ, origin, _rotation)
end

#Type promotion happens here
function TorusMantle(D,r_torus::TRTo, r_tube::TRTu, φ::Nothing, θ::TT, origin::PT, rotation::ROT) where {TRTo, TRTu, TT, PT, ROT}
    eltypes = _csg_get_promoted_eltype.((TRTo, TRTu, TT, PT, ROT))
    T = float(promote_type(eltypes...))
    TorusMantle{T}(D,r_torus, r_tube, φ, θ, origin, rotation)
end

function TorusMantle(D,r_torus::TRTo, r_tube::TRTu, φ::TP, θ::TT, origin::PT, rotation::ROT) where {TRTo, TRTu, TP, TT, PT, ROT}
    eltypes = _csg_get_promoted_eltype.((TRTo, TRTu, TP, TT, PT, ROT))
    T = float(promote_type(eltypes...))
    TorusMantle{T}(D,r_torus, r_tube, φ, θ, origin, rotation)
end

function TorusMantle(D::Symbol=:outwards;
    # define default parameters as Int to not influence type promotion
    r_torus = 1, 
    r_tube = 1,
    φ = nothing,
    θ = nothing,
    origin = zero(CartesianPoint{Int}), 
    rotation = one(SMatrix{3, 3, Int, 9})
) 
    TorusMantle(D,r_torus, r_tube, φ, θ, origin, rotation)
end

function TorusMantle{T}(D::Symbol=:outwards;
    r_torus = 1.0, 
    r_tube = 1.0,
    φ = nothing,
    θ= nothing,
    origin = zero(CartesianPoint{Float64}), 
    rotation = one(SMatrix{3, 3, Float64, 9})
) where {T}
    TorusMantle{T}(D,r_torus, r_tube, φ, θ, origin, rotation)
end

flip(t::TorusMantle{T,TP,TT,:inwards}) where {T,TP,TT} = 
TorusMantle{T,TP,TT,:outwards}(t.r_torus, t.r_tube, t.φ, t.θ, t.origin, t.rotation )
flip(t::TorusMantle{T,TP,TT,:outwards}) where {T,TP,TT} = 
TorusMantle{T,TP,TT,:inwards}(t.r_torus, t.r_tube, t.φ, t.θ, t.origin, t.rotation )

get_φ_limits(tm::TorusMantle{T,T}) where {T} = T(0), tm.φ
get_φ_limits(tm::TorusMantle{T,Nothing}) where {T} = T(0), T(2π)

get_θ_limits(tm::TorusMantle{T,<:Any,Tuple{T,T}}) where {T} = tm.θ[1], tm.θ[2]
get_θ_limits(tm::TorusMantle{T,<:Any,Nothing}) where {T} = T(0), T(2π)
get_θ_limits(tm::TorusMantle{T,<:Any,T}) where {T} = T(0), tm.θ

function normal(tm::TorusMantle{T,TP,TT,:outwards}, pt::CartesianPoint{T})::CartesianVector{T} where {T,TP,TT}
    pto = _transform_into_object_coordinate_system(pt, tm)
    cyl = CylindricalPoint(pto)
    ptt = CartesianPoint(CylindricalPoint{T}(tm.r_torus, cyl.φ, zero(T)))
    return CartesianVector(pt - _transform_into_global_coordinate_system(ptt, tm))
end
normal(tm::TorusMantle{T,TP,TT,:inwards}, pt::CartesianPoint{T}) where {T,TP,TT} = -normal(flip(tm), pt)

function vertices(tm::TorusMantle{T}, n_arc::Int)::Vector{CartesianPoint{T}} where {T}
    φMin, φMax = get_φ_limits(tm)
    θMin, θMax = get_θ_limits(tm)
    n_arcφ = _get_n_points_in_arc_φ(tm, n_arc) 
    n_arcθ = _get_n_points_in_arc_θ(tm, n_arc)
    
    θ = range(θMin, stop = θMax, length = n_arcθ + 1)
    φ = range(φMin, stop = φMax, length = n_arcφ + 1)
    
    [_transform_into_global_coordinate_system(CartesianPoint{T}((tm.r_torus + tm.r_tube*cos(θ))*cos(φ), (tm.r_torus + tm.r_tube*cos(θ))*sin(φ), tm.r_tube*sin(θ)), tm) for θ in θ for φ in φ]
end

function sample(tm::TorusMantle{T}, spacing::T)::Vector{CartesianPoint{T}} where {T}
    φMin, φMax = get_φ_limits(tm)
    θMin, θMax = get_θ_limits(tm)
    Δφ = abs(φMax - φMin)
    Δθ = abs(θMax - θMin)
    full2πφ = isnothing(tm.φ)
	full2πθ = isnothing(tm.θ)
	
	[_transform_into_global_coordinate_system(CartesianPoint{T}((tm.r_torus + tm.r_tube*cos(θ))*cos(φ), (tm.r_torus + tm.r_tube*cos(θ))*sin(φ), tm.r_tube*sin(θ)), tm) 
        for θ in range(θMin, stop = θMax - full2πθ*spacing/tm.r_tube, length = max(2,1+Int(ceil(Δθ*tm.r_tube/spacing)))) 
        for φ in range(φMin, stop = φMax - full2πφ*spacing/(tm.r_torus + tm.r_tube*cos(θ)), length = max(2,1+Int(ceil(Δφ*(tm.r_torus + tm.r_tube*cos(θ))/spacing))))]
end

function connections(tm::TorusMantle, n_arc::Int)::Vector{Vector{Int}}
    n_arcφ = _get_n_points_in_arc_φ(tm, n_arc) 
    n_arcθ = _get_n_points_in_arc_θ(tm, n_arc)
    [[i+(n_arcφ+1)*j,i+1+(n_arcφ+1)*j,i+1+(n_arcφ+1)*(j+1),i+(n_arcφ+1)*(j+1)] for j in 0:n_arcθ-1 for i in 1:n_arcφ]
end

function _get_hor_lines_idx(tm::TorusMantle{T}, n_arcθ::Int)::Vector{Int} where {T}
    θMin, θMax = get_θ_limits(tm) 
    if mod(θMax-θMin, T(2π)) == 0 
        l_0 = Int(floor(abs(θMin/(θMax-θMin))*n_arcθ)) + 1
        l_π = Int(floor(abs((θMin-π)/(θMax-θMin))*n_arcθ)) + 1
        [l_0,l_π]
    else
        [1, n_arcθ+1]
    end
end

function connections(tm::TorusMantle{T}, n_arc::Int, n_vert_lines::Int)::Vector{Vector{Int}} where {T}
    n_arcφ = _get_n_points_in_arc_φ(tm, n_arc) 
    n_arcθ = _get_n_points_in_arc_θ(tm, n_arc)
    vcircs = [[i*(n_arcφ+1)+c, (i+1)*(n_arcφ+1)+c] for c in _get_vert_lines_range(tm,n_arcφ,n_vert_lines) for i in 0:n_arcθ-1]
    hcircs =  [[i + (c-1)*(n_arcφ+1), i + (c-1)*(n_arcφ+1) + 1] for c in _get_hor_lines_idx(tm, n_arcθ) for i in 1:n_arcφ]
    append!(vcircs, hcircs)
end

get_label_name(::TorusMantle) = "Torus Mantle"

const FullTorusMantle{T,D} = TorusMantle{T,Nothing,Nothing,D}
const FullPhiTorusMantle{T,D} = TorusMantle{T,Nothing,Tuple{T,T},D}
const FullThetaTorusMantle{T,D} = TorusMantle{T,T,Nothing,D}

function lines(tm::TorusMantle{T,TP,Nothing}) where {T,TP}
    top_circ_origin = CartesianPoint{T}(zero(T), zero(T),  tm.r_tube)
    top_circ_origin = _transform_into_global_coordinate_system(top_circ_origin, tm)
    top_circ = Ellipse{T,T,TP}(tm.r_torus, tm.φ, top_circ_origin, tm.rotation)
    bot_circ_origin = CartesianPoint{T}(zero(T), zero(T), -tm.r_tube)
    bot_circ_origin = _transform_into_global_coordinate_system(bot_circ_origin, tm)
    bot_circ = Ellipse{T,T,TP}(tm.r_torus, tm.φ, bot_circ_origin, tm.rotation)
    inner_circ = Ellipse{T,T,TP}(tm.r_torus - tm.r_tube, tm.φ, tm.origin, tm.rotation)
    outer_circ = Ellipse{T,T,TP}(tm.r_torus + tm.r_tube, tm.φ, tm.origin, tm.rotation)
    φmin::T, φmax::T = isnothing(tm.φ) ? (0, π) : get_φ_limits(tm)
    tube_circ_1_origin = CartesianPoint(CylindricalPoint{T}(tm.r_torus, φmin, zero(T)))
    tube_circ_1_origin = _transform_into_global_coordinate_system(tube_circ_1_origin, tm)
    tube_circ_1 = Ellipse{T,T,Nothing}(tm.r_tube, tm.θ, tube_circ_1_origin, tm.rotation * RotZ(φmin) * RotX(T(π)/2))
    tube_circ_2_origin = CartesianPoint(CylindricalPoint{T}(tm.r_torus, φmax, zero(T)))
    tube_circ_2_origin = _transform_into_global_coordinate_system(tube_circ_2_origin, tm)
    tube_circ_2 = Ellipse{T,T,Nothing}(tm.r_tube, tm.θ, tube_circ_2_origin, tm.rotation * RotZ(φmax) * RotX(T(π)/2))
    (bot_circ, outer_circ, top_circ, inner_circ, tube_circ_1, tube_circ_2)
end 

function lines(tm::TorusMantle{T,TP,Tuple{T,T}}) where {T,TP}
    top_circ_origin = CartesianPoint{T}(zero(T), zero(T),  tm.r_tube)
    top_circ_origin = _transform_into_global_coordinate_system(top_circ_origin, tm)
    top_circ = Ellipse{T,T,TP}(tm.r_torus, tm.φ, top_circ_origin, tm.rotation)
    bot_circ_origin = CartesianPoint{T}(zero(T), zero(T), -tm.r_tube)
    bot_circ_origin = _transform_into_global_coordinate_system(bot_circ_origin, tm)
    bot_circ = Ellipse{T,T,TP}(tm.r_torus, tm.φ, bot_circ_origin, tm.rotation)
    inner_circ = Ellipse{T,T,TP}(tm.r_torus - tm.r_tube, tm.φ, tm.origin, tm.rotation)
    outer_circ = Ellipse{T,T,TP}(tm.r_torus + tm.r_tube, tm.φ, tm.origin, tm.rotation)
    φmin::T, φmax::T = isnothing(tm.φ) ? (0, π) : get_φ_limits(tm)
    tube_circ_1_origin = CartesianPoint(CylindricalPoint{T}(tm.r_torus, φmin, zero(T)))
    tube_circ_1_origin = _transform_into_global_coordinate_system(tube_circ_1_origin, tm)
    tube_circ_1 = Ellipse{T,T,T}(tm.r_tube, abs(-(tm.θ...)), tube_circ_1_origin, tm.rotation * RotZ(φmin) * RotX(T(π)/2) * RotZ(tm.θ[1]))
    tube_circ_2_origin = CartesianPoint(CylindricalPoint{T}(tm.r_torus, φmax, zero(T)))
    tube_circ_2_origin = _transform_into_global_coordinate_system(tube_circ_2_origin, tm)
    tube_circ_2 = Ellipse{T,T,T}(tm.r_tube, abs(-(tm.θ...)), tube_circ_2_origin, tm.rotation * RotZ(φmax) * RotX(T(π)/2) * RotZ(tm.θ[1]))
    (bot_circ, outer_circ, top_circ, inner_circ, tube_circ_1, tube_circ_2)
end 

extremum(tm::TorusMantle{T}) where {T} = tm.r_torus + tm.r_tube


# function _in(pt::CartesianPoint{T}, t::TorusMantle{T}; csgtol::T = csg_default_tol(T)) where {T}
#     return (isnothing(t.φ) || _in_angular_interval_closed(atan(pt.y, pt.x), t.φ, csgtol = csgtol)) &&
#            (isnothing(t.θ) || _in_angular_interval_closed(atan(pt.z, hypot(pt.x, pt.y) - t.r_torus), t.θ, csgtol = csgtol))
# end

"""
    intersection(tm::TorusMantle{T}, l::Line{T}) where {T}

Calculates the intersections of a `Line` with a `TorusMantle`.

## Arguments
* `cm::TorusMantle{T}`: The `TorusMantle`.
* `l::Line{T}`: The `Line`.

!!! note 
    The function will always return 4 CartesianPoint's.
    If the line just touches the mantle, the points will be the same. 
    If the line does not touch the mantle at all, the points will have NaN's as there coordinates.
"""
function intersection(tm::TorusMantle{T}, l::Line{T}) where {T}
    obj_l = _transform_into_object_coordinate_system(l, tm) # direction is not normalized
    
    let T2 = Float64 # will this work on 32bit-machines?

    L1::T2 = obj_l.origin.x |> T2
    L2::T2 = obj_l.origin.y |> T2
    L3::T2 = obj_l.origin.z |> T2
    D1::T2 = obj_l.direction.x |> T2
    D2::T2 = obj_l.direction.y |> T2
    D3::T2 = obj_l.direction.z |> T2

    R::T2 = tm.r_torus |> T2
    r::T2 = tm.r_tube  |> T2

    A::T2 = L1^2 + L2^2 + L3^2 + R^2 - r^2
    B::T2 = 2*(L1*D1 + L2*D2 + L3*D3)
    C::T2 = D1^2 + D2^2 + D3^2

    a::T2 = (2*B*C)
    b::T2 = (2*A*C + B^2 - 4*R^2*(D1^2 + D2^2))
    c::T2 = (2*A*B - 8*R^2*(L1*D1 + L2*D2))
    d::T2 = (A^2 - 4*R^2*(L1^2 + L2^2))

    # Solve: `solve (sqrt((L1 + λ*D1)^2 + (L2 + λ*D2)^2)-R)^2 + (L3 + λ*D3)^2 = r^2 for λ`

	# λ1, λ2, λ3, λ4 = roots_of_4th_order_polynomial(a, b, c, d) # That does not work for all combinations of a, b, c, d...
	# fallback to Polynomials.jl, which is slower... We should improve `roots_of_4th_order_polynomial`... 
    return broadcast(λ -> _transform_into_global_coordinate_system(obj_l.origin + λ * obj_l.direction, tm), 
        T.(real.(Polynomials.roots(Polynomial((d, c, b, a, C^2))))))
    # pts[.!in.(pts, Ref(tm))] .= Ref(CartesianPoint{T}(NaN,NaN,NaN)) # for Partial Tori: discard all intersection points that are not part of the Torus
    # return pts
    end
end

# These function compose cross sections of Tori at constant θ and ensure that the height is always positive
TorusThetaSurface(r::T, φ::TP, hZ::T, origin::CartesianPoint{T}, rotation::SMatrix{3,3,T,9}, ::Val{:flat}) where {T,TP} = EllipticalSurface{T,T,TP}(r,φ,origin,rotation)
TorusThetaSurface(r::T, φ::TP, hZ::T, origin::CartesianPoint{T}, rotation::SMatrix{3,3,T,9}, ::Val{:inwards}) where {T,TP} = ConeMantle{T,T,TP,:inwards}(r,φ,abs(hZ),origin,rotation)
TorusThetaSurface(r::T, φ::TP, hZ::T, origin::CartesianPoint{T}, rotation::SMatrix{3,3,T,9}, ::Val{:outwards}) where {T,TP} = ConeMantle{T,T,TP,:outwards}(r,φ,abs(hZ),origin,rotation)
TorusThetaSurface(r::Tuple{T,T}, φ::TP, hZ::T, origin::CartesianPoint{T}, rotation::SMatrix{3,3,T,9}, ::Val{:flat}) where {T,TP} = EllipticalSurface{T,Tuple{T,T},TP}(ordered(r...),φ,origin,rotation)
TorusThetaSurface(r::Tuple{T,T}, φ::TP, hZ::T, origin::CartesianPoint{T}, rotation::SMatrix{3,3,T,9}, ::Val{:inwards}) where {T,TP} = ConeMantle{T,Tuple{T,T},TP,:inwards}(hZ < 0 ? reverse(r) : r, φ, abs(hZ), origin, rotation)
TorusThetaSurface(r::Tuple{T,T}, φ::TP, hZ::T, origin::CartesianPoint{T}, rotation::SMatrix{3,3,T,9}, ::Val{:outwards}) where {T,TP} = ConeMantle{T,Tuple{T,T},TP,:outwards}(hZ < 0 ? reverse(r) : r, φ, abs(hZ), origin, rotation)


# """
#     roots_of_4th_order_polynomial(a::T, b::T, c::T, d::T, e::T)
# 
# Calculate the 4 (possible) roots of `x^4 + ax^3 + bx^2 + cx + d = 0`
# """
# function roots_of_4th_order_polynomial(a::T, b::T, c::T, d::T) where {T}
# 	#=
# 	using Polynomials
# 	A, a, b, c, d = 1.0, 0.2, -5.0, -0.1, 2.5
# 	p = Polynomial(SVector{5,T}([d, c, b, a, A]))
# 	xs = -2.5:0.01:2.5; plot(xs, map(x -> p(x), xs))
# 	roots(p)
# 
# 	rs = CSG.roots_of_4th_order_polynomial(a, b, c, d)
# 	vline!([rs...])
# 	@btime CSG.roots_of_4th_order_polynomial($a, $b, $c, $d)
# 	=#
# 
# 	# There are some issues for certain combinations of a, b, c, d, ...
# 
# 	term_1_1_1_1 = 2*b^3 - 9*a*b*c + 27*c^2 + 27*a^2*d - 72*b*d
# 	term_1_1_1_2 = b^2 - 3*a*c + 12*d
# 
# 	comp_term = -4 * term_1_1_1_2^3 + term_1_1_1_1^2
# 
# 	term_1_1_1_3 = term_1_1_1_1 + sqrt(Complex(comp_term))
# 	term_1_1_1 = 3*(term_1_1_1_3^(1/3))
# 
# 	term_1_1 = 2^(1/3) * term_1_1_1_2 / term_1_1_1
# 	term_1_2 = (term_1_1_1_3 / 54)^(1/3)
# 	term_1 = a^2/4 - 2b/3 + term_1_1 + term_1_2
# 
# 	if term_1 == 0
# 		term_6 = sqrt(Complex(b^2 - 4d))
# 		λ1 = sqrt(-term_6 - b)/sqrt(2)
# 		λ2 = -λ1
# 		λ3 = sqrt( term_6 - b)/sqrt(2)
# 		λ4 = -λ3
# 	else
# 		term_2_1 = (-a^3 + 4a*b - 8c) / (4*sqrt(term_1))
# 		term_2a = a^2/2 - 4b/3 - term_1_2 - term_1_2 - term_2_1
# 		term_2b = a^2/2 - 4b/3 - term_1_1 - term_1_2 + term_2_1
# 
# 		term3  = sqrt(term_1)/2
# 		term4a = sqrt(term_2a)/2 
# 		term4b = sqrt(term_2b)/2 
# 		term5 = -a/4	
# 
# 		λ1 = term5 - term3 - term4a 
# 		λ2 = term5 - term3 + term4a 
# 		λ3 = term5 + term3 - term4b 
# 		λ4 = term5 + term3 + term4b 
# 	end
# 	real(λ1), real(λ2), real(λ3), real(λ4)
# end


# function distance_to_surface(pt::AbstractCoordinatePoint{T}, t::TorusMantle{T, Nothing})::T where {T}
#     pcy = CylindricalPoint(pt)
#     return distance_to_line(PlanarPoint{T}(pcy.r,pcy.z), Arc(t))
# end

# function distance_to_surface(pt::AbstractCoordinatePoint{T}, t::TorusMantle{T, <:AbstractInterval})::T where {T}
#     pcy = CylindricalPoint(pt)
#     if _in_φ(pt, t.φ)
#         return distance_to_line(PlanarPoint{T}(pcy.r,pcy.z), Arc(t))
#     else
#         φMin::T, φMax::T, _ = get_φ_limits(t)
#         Δφ = pcy.φ - _φNear(pcy.φ, φMin, φMax)
#         d, r_on_plane = pcy.r .* sincos(Δφ)
#         return hypot(d, distance_to_line(PlanarPoint{T}(r_on_plane, pcy.z), Arc(t)))
#     end
# end


function _Δθ(θ1::T, θ2::T)::T where {T}
    δθ = mod(θ2 - θ1, T(2π))
    min(δθ, T(2π) - δθ)
end
_θNear(θ::Real, θMin::T, θMax::T) where {T} = _Δθ(T(θ),θMin) ≤ _Δθ(T(θ),θMax) ? θMin : θMax

function distance_to_surface(pt::AbstractCoordinatePoint{T}, t::TorusMantle{T, Nothing})::T where {T}

    pt_cyl  = CylindricalPoint(_transform_into_object_coordinate_system(CartesianPoint(pt), t))

    _in_φ(pt_cyl, t.φ) || throw(ArgumentError("φ-Partial Torus not implemented yet"))

    pt_r = pt_cyl.r - t.r_torus
    pt_z = pt_cyl.z
    pt_θ = atan(pt_z, pt_r)

    θMin::T, θMax::T = get_θ_limits(t)
    return if _in_angular_interval_closed(pt_θ, t.θ)
        abs(hypot(pt_r, pt_z) - t.r_tube)
    else
        sθNear, cθNear = sincos(_θNear(pt_θ, θMin, θMax))
        hypot(pt_r - t.r_tube * cθNear, pt_z - t.r_tube * sθNear)
    end
end