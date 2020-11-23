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

function Tube(;rMin::T = 0, rMax::T = 1, φMin::T = 0, φMax::T = 2π, zMin::T = -1/2, zMax::T = 1/2) where {T}
    Tube{T, rMin == 0 ? :full : :open, mod(φMax - φMin, 2π) == 0 ? :full : :open}(
        rMin, rMax, φMin, φMax, zMin, zMax
    )
end

in(p::CartesianPoint, t::Tube{<:Any, :full, :full}) = 
    (t.zMin <= p.z <= t.zMax) && ((p.x^2 + p.y^2) <= t.rMax^2)

in(p::CartesianPoint, t::Tube{<:Any, :open, :full}) = 
    (t.zMin <= p.z <= t.zMax) && (t.rMin^2 <= (p.x^2 + p.y^2) <= t.rMax^2)

in(p::CartesianPoint, t::Tube{<:Any, :full, :open}) = 
    (t.zMin <= p.z <= t.zMax) && ((p.x^2 + p.y^2) <= t.rMax^2) && (t.φMax <= mod(atan(pt.y, pt.x), 2π) <= t.φMax)

in(p::CartesianPoint, t::Tube{<:Any, :open, :open}) = 
    (t.zMin <= p.z <= t.zMax) && ((t.rMin^2 <= p.x^2 + p.y^2) <= t.rMax^2) && (t.φMax <= mod(atan(pt.y, pt.x), 2π) <= t.φMax)


