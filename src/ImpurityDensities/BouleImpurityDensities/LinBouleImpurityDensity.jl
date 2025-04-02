"""
    struct LinBouleImpurityDensity{T <: SSDFloat} <: AbstractImpurityDensity{T}

Impurity density model which assumes a linear gradient in impurity density in z.
 
## Fields
* `a:T`: impurity density values at the boule origin.
* `b::T`: linear slope in `z` direction.
* `det_z0::T`: z coordinate of the detector origin in boule coordinates. The z-dirrection of the detector is opposite to the z-direction of the boule coordinates.
"""

struct LinBouleImpurityDensity{T <: SSDFloat} <: AbstractImpurityDensity{T}
    # a + b*z
    a::T
    b::T 
    det_z0::T
end

function get_impurity_density(idm::LinBouleImpurityDensity, pt::AbstractCoordinatePoint{T})::T where {T}
    cpt = CartesianPoint(pt)
    z = cpt[3]
    idm.a .+ idm.b * (idm.det_z0 .- z)
end