struct RegularPolygon{N,T,TR} <: AbstractSurfacePrimitive{T}
    r::TR
    z::T
    function RegularPolygon( N::Integer,
                   ::Type{T},
                   r::Union{T, <:AbstractInterval{T}},
                   z::T) where {T}
        @assert typeof(N) <: Integer && N >= 3 "Failed constructing RegularPolygon with $(N) points in the polygon base (needs at least 3 points)."
        new{N,T,typeof(r)}(r, z)
    end
end

# Convenience functions
const RegularTriangle{T,TR} = RegularPolygon{3,T,TR}
const Square{T,TR}  = RegularPolygon{4,T,TR}
const RegularPentagon{T,TR} = RegularPolygon{5,T,TR}
const RegularHexagon{T,TR}  = RegularPolygon{6,T,TR}

RegularTriangle(args...) = RegularPolygon(3, args...)
Square(args...)  = RegularPolygon(4, args...)
RegularPentagon(args...) = RegularPolygon(5, args...)
RegularHexagon(args...)  = RegularPolygon(6, args...)

#Constructors
RegularPolygon(p::RegularPrism{N,T}; z = 0) where {N,T} = RegularPolygon(N, T, p.r, T(z))

function RegularPolygon(N::Integer; rInner = 0, rOuter = 1, z = 0)
    T = float(promote_type(typeof.((rInner, rOuter, z))...))
    r = rInner == 0 ? T(rOuter) : T(rInner)..T(rOuter)
    RegularPolygon(N, T, r, T(z))
end

RegularPolygon(N::Integer, rInner, rOuter, z) = RegularPolygon(N; rInner = rInner, rOuter = rOuter, z = z)

function RegularPolygon(N::Integer, rOuter::R, z::Z) where {R<:Real, Z<:Real}
    T = float(promote_type(R,Z))
    RegularPolygon( N, T, T(rOuter), T(z))
end

@inline in(p::CylindricalPoint, rp::RegularPolygon{N, T, <:Real}) where {N,T} = begin
    _isapprox_z(p, rp.z) && p.r * cos(T(π/N) - mod(p.φ, T(2π/N))) / cos(T(π/N)) <= rp.r
end

@inline in(p::CylindricalPoint, rp::RegularPolygon{N, T, <:AbstractInterval}) where {N,T} = begin
    _isapprox_z(p, rp.z) && p.r * cos(T(π/N) - mod(p.φ, T(2π/N))) / cos(T(π/N)) in rp.r
end

@inline in(p::CartesianPoint, rp::RegularPolygon) = in(CylindricalPoint(p), rp)

# Special case: CartesianPoint in RegularHexagon: use analytical formulas
@inline in(p::CartesianPoint, hp::RegularHexagon{T, <:Real}) where {T} =
    _isapprox_z(p, hp.z) && _in_y(p, hp.r * sqrt(T(3))/2) && _in_x(p, hp.r - abs(p.y)/sqrt(T(3)))

@inline in(p::CartesianPoint, hp::RegularHexagon{T, <:AbstractInterval{T}}) where {T} =
    _isapprox_z(p, hp.z) && abs(p.y) <= hp.r.right * sqrt(T(3))/2 &&
    (
        abs(p.y) >= hp.r.left * sqrt(T(3))/2 && _in_x(p, hp.r.right * T(0.5)) ||
        abs(p.x) in (hp.r.left - abs(p.y) /sqrt(T(3)))..(hp.r.right - abs(p.y)/sqrt(T(3)))
    )

get_r_limits(rp::RegularPolygon) = (_left_radial_interval(rp.r), _right_radial_interval(rp.r))
get_r_at_φ(rp::RegularPolygon{N,T}, φ::T) where {N,T} = get_r_limits(rp) .* T(cos(π/N)/cos(π/N - mod(φ, 2π/N)))

function get_y_at_x(r::T, corners::SVector{N, Tuple{T,T}}, g::CartesianTicksTuple{T}, x::T)::T where {N,T}
    for i in eachindex(corners)[1:end-1]
        s1,c1 = r .* corners[i]
        s2,c2 = r .* corners[i+1]
        if isapprox(c1, c2, atol = geom_atol_zero(T))
            return isapprox(x, c1, atol = geom_atol_zero(T)) ? max(s1,s2) : T(0)
        end
        t = (x-c1)/(c2-c1)
        if t in 0..1 return s1 + t * (s2-s1) end
    end
    NaN # If there is no point of RegularPolygon at x (used for the inner cutout of polygons)
end

function get_y_ticks(rMin::T, rMax::T, corners::SVector{N, Tuple{T,T}}, g::CartesianTicksTuple{T}, x::T) where {N,T}
    yMin::T = rMin == 0 ? NaN : get_y_at_x(rMin, corners, g, x)
    yMax::T = get_y_at_x(rMax, corners, g, x)
    isnan(yMin) ? _get_ticks(g.y, -yMax, yMax) : vcat(_get_ticks(g.y, -yMax, -yMin), _get_ticks(g.y, yMin, yMax))
end
    

function get_missing_x_at_y(r::T, corners::SVector{N, Tuple{T,T}}, g::CartesianTicksTuple{T}, y::T) where {N,T}
    x = []
    for i in eachindex(corners)[1:end-1]
        s1,c1 = r .* corners[i]
        s2,c2 = r .* corners[i+1]
        if ! isapprox(s1, s2, atol = geom_atol_zero(T))
            t = (y-s1)/(s2-s1)
            if t in 0..1 push!(x, c1 + t * (c2-c1)) end
        end
    end
    x
end


function sample(rp::RegularPolygon{N,T}, Nsamps::NTuple{3,Int} = (2,N+1,2))::Vector{CylindricalPoint{T}} where {N,T}
    rMin::T, rMax::T = get_r_limits(rp)
    samples = [
        CylindricalPoint{T}(r*cos(π/N)/cos(π/N - mod(φ, 2π/N)),φ,rp.z)
        for φ in (Nsamps[2] ≤ 1 ? 0 : range(0, 2π, length = Nsamps[2])) #ignore the N given in Nsamps and sample with N
        for r in (Nsamps[1] ≤ 1 ? rMax : range(rMin, rMax, length = Nsamps[1]))
        ]
end


function sample(rp::RegularPolygon{N,T}, g::CylindricalTicksTuple{T})::Vector{CylindricalPoint{T}} where {N,T}
    samples = [
        CylindricalPoint{T}(r,φ,rp.z)
        for φ in _get_ticks(g.φ, T(0), T(2π))
        for r in _get_ticks(g.r, get_r_at_φ(rp,φ)...)
    ]
end

function sample(rp::RegularPolygon{N,T}, g::CartesianTicksTuple{T})::Vector{CartesianPoint{T}} where {N,T}
    corners = SVector{N+1,Tuple{T,T}}(sincos.(T(2π/N) .* (0:N)))
    rMin::T, rMax::T = get_r_limits(rp)
    samples = vcat(
    [   # sample all points on the x-grid lines
        CartesianPoint{T}(x,y,rp.z)
        for x in _get_ticks(g.x, rMax * minimum(broadcast(p -> p[2], corners)), rMax * maximum(broadcast(p -> p[2], corners)))
        for y in get_y_ticks(rMin, rMax, corners, g, x)
    ] , 
    [   # sample the missing points on y-grid lines
        CartesianPoint{T}(x,y,rp.z)
        for y in _get_ticks(g.y, rMax * minimum(broadcast(p -> p[1], corners)), rMax * maximum(broadcast(p -> p[1], corners)))
        for x in vcat(get_missing_x_at_y(rMin, corners, g, y), get_missing_x_at_y(rMax, corners, g, y))
    ]
    )
end

function distance_to_surface(point::AbstractCoordinatePoint{T}, rp::RegularPolygon{N,T})::T where {N,T}
    pcy = CylindricalPoint(point)
    rMin::T, rMax::T = get_r_limits(rp)
    φ = mod(pcy.φ, T(2π/N))
    r_poly = pcy.r * cos(T(π/N) - φ) / cos(T(π/N))
    if rMin ≤ r_poly ≤ rMax
        return abs(pcy.z - rp.z)
    else
        sn, cn = sincos(2π/N)
        r = r_poly < rMin ? rMin : rMax
        line = LineSegment(T, PlanarPoint{T}(r, 0), PlanarPoint{T}(r*cn, r*sn))
        sφ, cφ = sincos(φ)
        d = distance_to_line(PlanarPoint{T}(pcy.r*cφ, pcy.r*sφ), line)
        return hypot(d, pcy.z - rp.z)    
    end
end