struct Cone{T,TR,TP,TZ} <: AbstractVolumePrimitive{T}
    rtop::TR
    rbot::TR
    φ::TP
    z::TZ
    
    function Cone( ::Type{T},
                   rtop::Union{T, <:AbstractInterval{T}},
                   rbot::Union{T, <:AbstractInterval{T}},
                   φ::Union{Nothing, <:AbstractInterval{T}},
                   z::Union{T, <:AbstractInterval{T}}) where {T}
        @assert typeof(rtop) == typeof(rbot) "r-Intervals in Cone are not type stable"
        new{T,typeof(rtop),typeof(φ),typeof(z)}(rtop, rbot, φ, z)
    end
end

#Constructors
function Cone(;rtopMin = 0, rtopMax = 1, rbotMin = 0, rbotMax = 1, φMin = 0, φMax = 2π, zMin = -1/2, zMax = 1/2)
    T = promote_type(typeof.((rtopMin, rtopMax, rbotMin, rbotMax, φMin, φMax, zMin, zMax))...)
    rtop, rbot = if rtopMin == rbotMin == 0 
        T(rtopMax), T(rbotMax)
    else
        T(rtopMin)..T(rtopMax), T(rbotMin)..T(rbotMax)
    end
    φ = mod(φMax - φMin, 2π) == 0 ? nothing : T(φMin)..T(φMax)
    z = zMax == -zMin ? T(zMax) : T(zMin)..T(zMax)
    rtop == rbot ? Tube( T, rtop, φ, z) : Cone( T, rtop, rbot, φ, z)
end
Cone(rtopMin, rtopMax, rbotMin, rbotMax, φMin, φMax, zMin, zMax) = Cone(;rtopMin, rtopMax, rbotMin, rbotMax, φMin, φMax, zMin, zMax)

function Cone(rtop::R1, rbot::R2, height::H) where {R1, R2, H}
    T = promote_type(R1, R2, H)
    Cone( T, T(rtop), T(rbot), nothing, T(height/2))
end
    
function get_r_at_z(c::Cone{T, <:AbstractInterval, <:Any, <:Any}, z::Real) where {T}
    r1::T = get_r_at_z(c.rtop.left, c.rbot.left, c.z, z)
    r2::T = get_r_at_z(c.rtop.right, c.rbot.right, c.z, z)
    r1..r2
end

function get_r_at_z(c::Cone{T, <:Real, <:Any, <:Any}, z::Real) where {T}    
    get_r_at_z(c.rtop, c.rbot, c.z, z)
end

@inline get_r_at_z(rtop::TR, rbot::TR, cz::TZ, z::Real) where {TR, TZ} = (rtop - rbot) * (z - _left(cz)) / _width(cz) + rbot

in(p::AbstractCoordinatePoint, c::Cone{<:Any, <:Any, Nothing, <:Any}) =
    _in_z(p, c.z) && _in_r(p, get_r_at_z(c, p.z))


in(p::AbstractCoordinatePoint, c::Cone{<:Any, <:Any, <:AbstractInterval, <:Any}) =
    _in_z(p, c.z) && _in_φ(p, c.φ) && _in_r(p, get_r_at_z(c, p.z))