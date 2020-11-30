struct HexagonalPrism{T,TR,TZ} <: AbstractVolumePrimitive{T}
    r::TR
    z::TZ
    
    function HexagonalPrism( ::Type{T},
                             r::Union{T, <:AbstractInterval{T}},
                             z::Union{T, <:AbstractInterval{T}}) where {T}
        new{T,typeof(r),typeof(z)}(r,z)
    end
end

# Constructors
function HexagonalPrism(; rInner = 0, rOuter = 1, zMin = -1, zMax = 1)
    T = float(promote_type(typeof.((rInner, rOuter, zMin, zMax))...))
    r = rInner == 0 ? T(rOuter) : T(rInner)..T(rOuter)
    z = zMax == -zMin ? T(zMax) : T(zMin)..T(zMax)
    HexagonalPrism( T, r, z)
end
HexagonalPrism(rInner, rOuter, zMin, zMax) = HexagonalPrism(; rInner, rOuter, zMin, zMax)

function HexagonalPrism(rOuter::R, height::H) where {R<:Real, H<:Real}
    T = float(promote_type(R,H))
    HexagonalPrism( T, T(rOuter), T(height/2))
end

in(p::CartesianPoint, hp::HexagonalPrism{T, <:Real, <:Any}) where {T} =
    _in_z(p, hp.z) && _in_x(p, hp.r * sqrt(T(3))/2) && _in_y(p, hp.r - abs(p.x)/sqrt(T(3)))
    
in(p::CartesianPoint, hp::HexagonalPrism{T, <:AbstractInterval{T}, <:Any}) where {T} = 
    _in_z(p, hp.z) && abs(p.x) <= hp.r.right * sqrt(T(3))/2 && 
    (
        abs(p.x) >= hp.r.left * sqrt(T(3))/2 && _in_y(p, hp.r.right * T(0.5)) || 
        abs(p.y) in (hp.r.left - abs(p.x) /sqrt(T(3)))..(hp.r.right - abs(p.x)/sqrt(T(3)))
    )

in(p::CylindricalPoint, hp::HexagonalPrism{T, <:Real, <:Any}) where {T} = begin
    φr::T = mod(p.φ, T(π/3))
    φ0::T = min(φr, T(π/3) - φr)
    _in_z(p, hp.z) && abs(p.r * T(2)/sqrt(T(3)) * cos(φ0)) <= hp.r
end

in(p::CylindricalPoint, hp::HexagonalPrism{T, <:AbstractInterval, <:Any}) where {T} = begin
    φr::T = mod(p.φ, T(π)/3)
    φ0::T = min(φr, T(π)/3 - φr)
    _in_z(p, hp.z) && p.r * T(2)/sqrt(T(3)) * cos(φ0) in hp.r
end

#=
#The important points will change with transformations. This will be reworked

# R
function get_important_points(hp::HexagonalPrism{T, <:Any, <:Any}, ::Val{:r})::Vector{T} where {T}
    return geom_round.(T[-_right_linear_interval_linear_interval(hp.r), -sqrt(T(3))/2*_right_linear_interval(hp.r), -_left_linear_interval(hp.r), -sqrt(T(3))/2*_left_linear_interval(hp.r), 0, sqrt(T(3))/2*_left_linear_interval(hp.r), _left_linear_interval(hp.r), sqrt(T(3))/2*_right_linear_interval(hp.r), _right_linear_interval(hp.r)])
end

# PHI
function get_important_points(hp::HexagonalPrism{T, <:Any, <:Any}, ::Val{:φ})::Vector{T} where {T}
    return geom_round.(T[π/6, 3π/6, 5π/6, 7π/6, 9π/6, 11π/6])
end

# Z
function get_important_points(hp::HexagonalPrism{T, <:Any, <:Any}, ::Val{:z})::Vector{T} where {T}
    return geom_round.(T[_left_linear_interval(hp.z), _right_linear_interval(hp.z)])
end

# X
function get_important_points(hp::HexagonalPrism{T, <:Any, <:Any}, ::Val{:x})::Vector{T} where {T}
    return geom_round.(T[-sqrt(T(3))/2*_right_linear_interval(hp.r), -sqrt(T(3))/2*_left_linear_interval(hp.r), 0, sqrt(T(3))/2*_left_linear_interval(hp.r), sqrt(T(3))/2*_right_linear_interval(hp.r)])
end

# Y
function get_important_points(hp::HexagonalPrism{T, <:Any, <:Any}, ::Val{:y})::Vector{T} where {T}
    return geom_round.(T[-_right_linear_interval(hp.r), -_right_linear_interval(hp.r)/2, -_left_linear_interval(hp.r), -_left_linear_interval(hp.r)/2, 0, _left_linear_interval(hp.r)/2, _left_linear_interval(hp.r), _right_linear_interval(hp.r)/2, _right_linear_interval(hp.r)])
end
=#
    

function sample(hp::HexagonalPrism{T, <:Any, <:Any}, stepsize::Vector{T}) where {T}
    samples = CartesianPoint{T}[]
    for x in -_right_linear_interval(hp.r) : stepsize[1] : _right_linear_interval(hp.r)
        for y in -_right_linear_interval(hp.r) : stepsize[2] : _right_linear_interval(hp.r)
            if CartesianPoint{T}(x, y, 0) in hp
                for z in _left_linear_interval(hp.z) : stepsize[3] : _right_linear_interval(hp.z)
                    push!(samples, CartesianPoint{T}(x, y, z))
                end
            end
        end
    end
    return samples
end
