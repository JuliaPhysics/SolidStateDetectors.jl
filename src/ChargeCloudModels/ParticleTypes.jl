"""
    abstract type ParticleType
        
Type of a particle that deposits energy in a [`Semiconductor`](@ref).
Currently defined are `Alpha`, `Beta` and `Gamma`.

`ParticleType` is used to determine the radius of an [`NBodyChargeCloud`](@ref),
where the default radius for `Alpha` is 0.1mm and the default radius for `Beta`
and `Gamma` is 0.5mm.
"""
abstract type ParticleType end
abstract type Alpha <: ParticleType end
abstract type Beta <: ParticleType end
abstract type Gamma <: ParticleType end

radius_guess(charge::T, ::Type{Alpha}) where {T} = T(0.0001)
radius_guess(charge::T, ::Type{Beta}) where {T} = T(0.0005)
radius_guess(charge::T, ::Type{Gamma}) where {T} = T(0.0005) 


struct SSDSource
  particle_type::String
  position::SolidStateDetectors.CartesianPoint
  direction::Union{Symbol, SolidStateDetectors.CartesianVector}
end

export SSDSource