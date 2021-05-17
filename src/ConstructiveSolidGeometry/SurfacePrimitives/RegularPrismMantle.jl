struct RegularPrismMantle{N,T,TZ} <: AbstractSurfacePrimitive{T}
    r::T
    z::TZ
    function RegularPrismMantle( N::Integer,
                   ::Type{T},
                   r::T,
                   z::Union{T, <:AbstractInterval{T}}) where {T}
        @assert typeof(N) <: Integer && N >= 3 "Failed constructing RegularPrismMantle with $(N) points in the polygon base (needs at least 3 points)."
        new{N,T,typeof(z)}(r, z)
    end
end

# Convenience functions
const TriangularPrismMantle{T,TZ} = RegularPrismMantle{3,T,TZ}
const SquarePrismMantle{T,TZ}  = RegularPrismMantle{4,T,TZ}
const PentagonalPrismMantle{T,TZ} = RegularPrismMantle{5,T,TZ}
const HexagonalPrismMantle{T,TZ}  = RegularPrismMantle{6,T,TZ}

TriangularPrismMantle(args...) = RegularPrismMantle(3, args...)
SquarePrismMantle(args...)  = RegularPrismMantle(4, args...)
PentagonalPrismMantle(args...) = RegularPrismMantle(5, args...)
HexagonalPrismMantle(args...)  = RegularPrismMantle(6, args...)

#Constructors
RegularPrismMantle(p::RegularPrism{N,T}; r = 1) where {N,T} = RegularPrismMantle(N, T, T(r), p.z)

function RegularPrismMantle(N::Integer; r = 1, zMin = -1, zMax = 1)
    T = float(promote_type(typeof.((r, zMin, zMax))...))
    z = zMax == -zMin ? T(zMax) : T(zMin)..T(zMax)
    RegularPrismMantle(N, T, T(r), z)
end
RegularPrismMantle(N::Integer, r, zMin, zMax) = RegularPrismMantle(N; r = r, zMin = zMin, zMax = zMax)

function RegularPrismMantle(N::Integer, r::R, height::H) where {R<:Real, H<:Real}
    T = float(promote_type(R,H))
    RegularPrismMantle( N, T, T(r), T(height/2))
end


@inline in(p::CylindricalPoint, rp::RegularPrismMantle{N, T}) where {N,T} = begin
    _in_z(p, rp.z) && isapprox(p.r * cos(T(π/N) - mod(p.φ, T(2π/N))) / cos(T(π/N)), rp.r, atol = geom_atol_zero(T))
end

@inline in(p::CartesianPoint, rp::RegularPrismMantle) = in(CylindricalPoint(p), rp)

@inline in(p::CartesianPoint, hp::HexagonalPrismMantle{T}) where {T} = begin
    tol = geom_atol_zero(T)
    _in_z(p, hp.z) && abs(p.y) ≤ hp.r * sqrt(T(3))/2 &&
    (
        ( isapprox(abs(p.y), hp.r * sqrt(T(3))/2, atol = tol) && _in_x(p, hp.r * T(0.5)) ) ||
        ( isapprox(abs(p.x), hp.r - abs(p.y)/sqrt(T(3)), atol = tol) &&  T(0.5) ≤ abs(p.x) ≤ hp.r )
    )
end

get_z_limits(rp::RegularPrismMantle) = (_left_linear_interval(rp.z), _right_linear_interval(rp.z))
get_r_at_φ(rp::RegularPrismMantle{N,T}, φ::T) where {N,T} = rp.r * T(cos(π/N)/cos(π/N - mod(φ, 2π/N)))

function sample(rp::RegularPrismMantle{N,T}, Nsamps::NTuple{3,Int} = (2,N+1,2))::Vector{CylindricalPoint{T}} where {N,T}
    zMin::T, zMax::T = get_z_limits(rp)
    samples = [
        CylindricalPoint{T}(rp.r*cos(π/N)/cos(π/N - mod(φ, 2π/N)),φ,z)
        for z in (Nsamps[3] ≤ 1 ? zMin : range(zMin, zMax, length = Nsamps[3]))
        for φ in (Nsamps[2] ≤ 1 ? 0 : range(0, 2π, length = Nsamps[2]))
        ]
end


function sample(rp::RegularPrismMantle{N,T}, g::CylindricalTicksTuple{T})::Vector{CylindricalPoint{T}} where {N,T}
    samples = [
        CylindricalPoint{T}(get_r_at_φ(rp,φ),φ,z)
        for φ in _get_ticks(g.φ, T(0), T(2π))
        for z in get_z_ticks(rp, g)
    ]
end

function sample(rp::RegularPrismMantle{N,T}, g::CartesianTicksTuple{T})::Vector{CartesianPoint{T}} where {N,T}
    corners = SVector{N+1,Tuple{T,T}}(sincos.(T(2π/N) .* (0:N)))
    samples = vcat(
    [   # sample all points on the x-grid lines
        CartesianPoint{T}(x,y,z)
        for x in _get_ticks(g.x, rp.r * minimum(broadcast(p -> p[2], corners)), rp.r * maximum(broadcast(p -> p[2], corners)))
        for y in (-get_y_at_x(rp.r, corners, g, x), get_y_at_x(rp.r, corners, g, x))
        for z in get_z_ticks(rp, g)
    ] , 
    [   # sample the missing points on y-grid lines
        CartesianPoint{T}(x,y,z)
        for y in _get_ticks(g.y, rp.r * minimum(broadcast(p -> p[1], corners)), rp.r * maximum(broadcast(p -> p[1], corners)))
        for x in get_missing_x_at_y(rp.r, corners, g, y)
        for z in get_z_ticks(rp, g)
    ]
    )
end

