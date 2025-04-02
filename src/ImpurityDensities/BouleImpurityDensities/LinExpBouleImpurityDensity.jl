"""
    struct LinExpBouleImpurityDensity{T <: SSDFloat} <: AbstractImpurityDensity{T}

Impurity density model which assumes a linear + exponential component in impurity density in z.
 
## Fields
* `a:T`: impurity density values at the boule origin.
* `b::T`: linear slope in `z` direction.
* `c::T`: exponential coefficient.
* `L::T`: exponential offset.
* `tau::T`: exponential scale factor.
* `det_z0::T`: z coordinate of the detector origin in boule coordinates. The z-dirrection of the detector is opposite to the z-direction of the boule coordinates.
"""

struct LinExpBouleImpurityDensity{T <: SSDFloat} <: AbstractImpurityDensity{T}
    # a + b*z + c*exp((z-L)/tau) -> needs at least 4 points
    a::T
    b::T 
    c::T 
    L::T
    tau::T
    det_z0::T
end

function get_impurity_density(idm::LinExpBouleImpurityDensity, pt::AbstractCoordinatePoint{T})::T where {T}
    cpt = CartesianPoint(pt)
    z = cpt[3]
    idm.a .+ idm.b * (idm.det_z0 .- z) .+ idm.c * exp.((idm.det_z0 .- z .- idm.L)/idm.tau)
end