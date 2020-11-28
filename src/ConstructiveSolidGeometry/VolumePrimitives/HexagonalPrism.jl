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
    T = promote_type(typeof.((rInner, rOuter, zMin, zMax))...)
    r = rInner == 0 ? T(rOuter) : T(rInner)..T(rOuter)
    z = zMax == -zMin ? T(zMax) : T(zMin)..T(zMax)
    HexagonalPrism( T, r, z)
end
HexagonalPrism(rInner, rOuter, zMin, zMax) = HexagonalPrism(; rInner, rOuter, zMin, zMax)

function HexagonalPrism(rOuter::R, height::H) where {R<:Real, H<:Real}
    T = promote_type(R,H)
    HexagonalPrism( T, T(rOuter), T(height/2))
end

in(p::CartesianPoint, hp::HexagonalPrism{<:Any, <:Real, <:Any}) =
    _in_z(p, hp.z) && _in_x(p, hp.r * sqrt(3)/2) && _in_y(p, hp.r - abs(p.x) /sqrt(3))
    
in(p::CartesianPoint, hp::HexagonalPrism{<:Any, <:AbstractInterval{T}, <:Any}) where {T} = 
    _in_z(p, hp.z) && abs(p.x) <= hp.r.right * sqrt(3)/2 && 
    (
        abs(p.x) >= hp.r.left * sqrt(3)/2 && _in_y(p, hp.r.right * T(0.5)) || 
        abs(p.y) in (hp.r.left - abs(p.x) /sqrt(3))..(hp.r.right - abs(p.x)/sqrt(3))
    )

in(p::CylindricalPoint{T}, hp::HexagonalPrism{<:Any, <:Real, <:Any}) where {T} = begin
    φr::T = mod(p.φ, T(π/3))
    φ0::T = min(φr, T(π/3) - φr)
    _in_z(p, hp.z) && abs(p.r * T(2/sqrt(3)) * cos(φ0)) <= hp.r
end

in(p::CylindricalPoint{T}, hp::HexagonalPrism{<:Any, <:AbstractInterval, <:Any}) where {T} = begin
    φr::T = mod(p.φ, T(π/3))
    φ0::T = min(φr, T(π/3) - φr)
    _in_z(p, hp.z) && p.r * T(2/sqrt(3)) * cos(φ0) in hp.r
end
    
# R
function get_important_points(hp::HexagonalPrism{T, <:Any, <:Any}, ::Val{:r})::Vector{T} where {T}
    return geom_round.(T[-_right(hp.r), -sqrt(3)/2*_right(hp.r), -_left(hp.r), -sqrt(3)/2*_left(hp.r), 0., sqrt(3)/2*_left(hp.r), _left(hp.r), sqrt(3)/2*_right(hp.r), _right(hp.r)])
end

# PHI
function get_important_points(hp::HexagonalPrism{T, <:Any, <:Any}, ::Val{:φ})::Vector{T} where {T}
    return geom_round.(T[π/6, 3π/6, 5π/6, 7π/6, 9π/6, 11π/6])
end

# Z
function get_important_points(hp::HexagonalPrism{T, <:Any, <:Any}, ::Val{:z})::Vector{T} where {T}
    return geom_round.(T[_left(hp.z), _right(hp.z)])
end

# X
function get_important_points(hp::HexagonalPrism{T, <:Any, <:Any}, ::Val{:x})::Vector{T} where {T}
    return geom_round.(T[-sqrt(3)/2*_right(hp.r), -sqrt(3)/2*_left(hp.r), 0., sqrt(3)/2*_left(hp.r), sqrt(3)/2*_right(hp.r)])
end

# Y
function get_important_points(hp::HexagonalPrism{T, <:Any, <:Any}, ::Val{:y})::Vector{T} where {T}
    return geom_round.(T[-_right(hp.r), -_right(hp.r)/2, -_left(hp.r), -_left(hp.r)/2, 0., _left(hp.r)/2, _left(hp.r), _right(hp.r)/2, _right(hp.r)])
end
    

function sample(hp::HexagonalPrism{T, <:Any, <:Any}, stepsize::Vector{<:Real}) where {T}
    samples = CartesianPoint{T}[]
    for x in -_right(hp.r) : stepsize[1] : _right(hp.r)
        for y in -_right(hp.r) : stepsize[2] : _right(hp.r)
            if CartesianPoint{T}(x, y, 0) in hp
                for z in _left(hp.z) : stepsize[3] : _right(hp.z)
                    push!(samples, CartesianPoint{T}(x, y, z))
                end
            end
        end
    end
    return samples
end