struct ConalPlane{T,TR,TP,TZ} <: AbstractSurfacePrimitive{T}
    r::TR #if tupple trapezoid/triangle, or else rectangle
    φ::TP
    z::TZ
    #if r is a Tuple, the first entry refers to the r-interval at the bottom, the second one to the r-interval at the top
    function ConalPlane( ::Type{T},
                   r::Union{T, <:AbstractInterval{T}, Tuple{T,T}, Tuple{I,I}},
                   φ::T,
                   z::Union{T, <:AbstractInterval{T}}) where {T, I<:AbstractInterval{T}}
        new{T,typeof(r),T,typeof(z)}(r, φ, z)
    end
end

ConalPlane(c::Cone{T}; φ = 0) where {T} = ConalPlane( T, c.r, T(mod(φ,2π)), c.z)

function ConalPlane(;rbotMin = 0, rbotMax = 1, rtopMin = 0, rtopMax = 1, φ = 0, zMin = -1/2, zMax = 1/2)
    T = float(promote_type(typeof.((rtopMin, rtopMax, rbotMin, rbotMax, φ, zMin, zMax))...))
    c = Cone(rbotMin, rbotMax, rtopMin, rtopMax, 0, 2π, zMin, zMax)
    ConalPlane( T, c.r, T(φ), c.z)
end
ConalPlane(rbotMin, rbotMax, rtopMin, rtopMax, φ, zMin, zMax) = ConalPlane(; rbotMin = rbotMin, rbotMax = rbotMax, rtopMin = rtopMin, rtopMax = rtopMax, φ = φ, zMin = zMin, zMax = zMax)

get_r_at_z(c::ConalPlane{T}, z::Real) where {T} = get_r_at_z(Cone(T, c.r, nothing, c.z), z::Real)
get_r_limits(c::ConalPlane{T}) where {T} = get_r_limits(Cone(T, c.r, nothing, c.z))
get_z_limits(c::ConalPlane{T}) where {T} = (_left_linear_interval(c.z), _right_linear_interval(c.z))

in(p::AbstractCoordinatePoint, c::ConalPlane) =
    _in_z(p, c.z) && _eq_φ(p, c.φ) && _in_cyl_r(p, get_r_at_z(c, p.z))

#function sample(c::ConalPlane{T}, step::Quantity{<:Real, Unitful.𝐋}) where {T}
function sample(c::ConalPlane{T}, step::Real) where {T}
    zMin::T, zMax::T = get_z_limits(c)
    #step = T(ustrip(uconvert(u"m", step)))
    samples = [
        CylindricalPoint{T}(r,c.φ,z)
        for z in zMin:step:zMax
        for r in _left_radial_interval(get_r_at_z(c, z)):step:_right_radial_interval(get_r_at_z(c, z))
    ]
end

function sample(c::ConalPlane{T}, Nsamps::NTuple{3,Int}) where {T}
    zMin::T, zMax::T = get_z_limits(c)
    samples = [
        CylindricalPoint{T}(r,c.φ,z)
        for z in (Nsamps[3] ≤ 1 ? zMin : range(zMin, zMax, length = Nsamps[3]))
        for r in (Nsamps[3] ≤ 1 ? _left_radial_interval(get_r_at_z(c, z)) : range(_left_radial_interval(get_r_at_z(c, z)), _right_radial_interval(get_r_at_z(c, z)), length = Nsamps[1]))
    ]
end

function get_vertices(c::ConalPlane{T}) where {T}
    rbotMin::T, rbotMax::T, rtopMin::T, rtopMax::T = get_r_limits(c)
    zMin::T, zMax::T = get_z_limits(c)
    [CartesianPoint{T}(rbotMin * cos(c.φ), rbotMin * sin(c.φ), zMin),
    CartesianPoint{T}(rbotMax * cos(c.φ), rbotMax * sin(c.φ), zMin),
    CartesianPoint{T}(rtopMax * cos(c.φ), rtopMax * sin(c.φ), zMax),
    CartesianPoint{T}(rtopMin * cos(c.φ), rtopMin * sin(c.φ), zMax)]
end
