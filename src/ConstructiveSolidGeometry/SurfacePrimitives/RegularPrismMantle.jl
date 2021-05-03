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
const TriangularPrismMantle{T,TR} = RegularPrismMantle{3,T,TR}
const SquarePrismMantle{T,TR}  = RegularPrismMantle{4,T,TR}
const PentagonalPrismMantle{T,TR} = RegularPrismMantle{5,T,TR}
const HexagonalPrismMantle{T,TR}  = RegularPrismMantle{6,T,TR}

TriangularPrismMantle(args...) = RegularPrismMantle(3, args...)
SquarePrismMantle(args...)  = RegularPrismMantle(4, args...)
PentagonalPrismMantle(args...) = RegularPrismMantle(5, args...)
HexagonalPrismMantle(args...)  = RegularPrismMantle(6, args...)

#Constructors
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
    _in_z(p, rp.z) && is_approx(p.r * cos(T(π/N) - mod(p.φ, T(2π/N))) / cos(T(π/N)), rp.r, atol = geom_atol_zero(T))
end

@inline in(p::CartesianPoint, rp::RegularPrismMantle) = in(CylindricalPoint(p), rp)

get_z_limits(rp::RegularPrismMantle) = (_left_linear_interval(rp.z), _right_linear_interval(rp.z))

function sample(rp::RegularPrismMantle{N,T}, Nsamps::NTuple{3,Int} = (2,N+1,2))::Vector{CylindricalPoint{T}} where {N,T}
    zMin::T, zMax::T = get_z_limits(rp)
    samples = [
        CylindricalPoint{T}(rp.r*cos(π/N)/cos(π/N - mod(φ, 2π/N)),φ,z)
        for z in (Nsamps[3] ≤ 1 ? zMin : range(zMin, zMax, length = Nsamps[3]))
        for φ in (Nsamps[2] ≤ 1 ? 0 : range(0, 2π, length = Nsamps[2]))
        ]
end
