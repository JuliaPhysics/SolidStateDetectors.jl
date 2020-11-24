"""
    struct Tube{T,R,PHI} <: AbstractVolumePrimitive{T}

`R` and `PHI` ∈ `[:full, :open]`
`R == :full` -> `rMin == 0`
`PHI == :full` -> `φMin == 0 && φMax == 2π`
"""
struct Tube{T,R,PHI} <: AbstractVolumePrimitive{T}
    rMin::T
    rMax::T
    φMin::T
    φMax::T 
    zMin::T
    zMax::T
end

# Constructors
function Tube(;rMin = 0, rMax = 1, φMin = 0, φMax = 2π, zMin = -1/2, zMax = 1/2) 
    T = promote_type(typeof.((rMin, rMax, φMin, φMax, zMin, zMax))...)
    Tube{T, rMin == 0 ? :full : :open, mod(φMax - φMin, 2π) == 0 ? :full : :open}(
        rMin, rMax, φMin, φMax, zMin, zMax
    )
end
Tube(rMin, rMax, φMin, φMax, zMin, zMax) = Tube(;rMin, rMax, φMin, φMax, zMin, zMax)    
Tube(radius::T, height::T) where {T} = Tube(rMax = radius, zMin = -height/2, zMax = height/2)

in(p::CartesianPoint, t::Tube{<:Any, :full, :full}) = 
    (t.zMin <= p.z <= t.zMax) && ((p.x^2 + p.y^2) <= t.rMax^2)

in(p::CartesianPoint, t::Tube{<:Any, :open, :full}) = 
    (t.zMin <= p.z <= t.zMax) && (t.rMin^2 <= (p.x^2 + p.y^2) <= t.rMax^2)

in(p::CartesianPoint, t::Tube{<:Any, :full, :open}) = 
    (t.zMin <= p.z <= t.zMax) && ((p.x^2 + p.y^2) <= t.rMax^2) && (t.φMax <= mod(atan(pt.y, pt.x), 2π) <= t.φMax)

in(p::CartesianPoint, t::Tube{<:Any, :open, :open}) = 
    (t.zMin <= p.z <= t.zMax) && ((t.rMin^2 <= p.x^2 + p.y^2) <= t.rMax^2) && (t.φMax <= mod(atan(pt.y, pt.x), 2π) <= t.φMax)


