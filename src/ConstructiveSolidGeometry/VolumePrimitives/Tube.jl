struct Tube{T,TR,TP,TZ} <: AbstractVolumePrimitive{T}
    r::TR
    φ::TP
    z::TZ

    function Tube(  ::Type{T},
                    r::Union{T, <:AbstractInterval{T}}, 
                    φ::Union{Nothing, <:AbstractInterval{T}}, 
                    z::Union{T, <:AbstractInterval{T}}) where {T}
        new{T,typeof(r),typeof(φ),typeof(z)}(r, φ, z)
    end
end

# Constructors
function Tube(;rMin = 0, rMax = 1, φMin = 0, φMax = 2π, zMin = -1/2, zMax = 1/2) 
    T = promote_type(typeof.((rMin, rMax, φMin, φMax, zMin, zMax))...)
    φ = mod(φMax - φMin, 2π) == 0 ? nothing : T(φMin)..T(φMax)
    z = zMax == -zMin ? T(zMax) : T(zMin)..T(zMax)
    Tube( T, T(rMin)..T(rMax), φ, z)
end
Tube(rMin, rMax, φMin, φMax, zMin, zMax) = Tube(;rMin, rMax, φMin, φMax, zMin, zMax)    

function Tube(radius::R, height::H) where {R,H}
    T = promote_type(R,H)
    Tube( T, T(radius), nothing, T(height/2))
end

in(p::CartesianPoint, t::Tube{<:Any, <:Any, Nothing, <:Any}) = 
    _in_z(p, t.z) && _in_r(p, t.r)

in(p::CartesianPoint, t::Tube{<:Any, <:Any, <:AbstractInterval, <:Any}) = 
    _in_z(p, t.z) && _in_r(p, t.r) && _in_φ(p, t.φ)
