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

"""
    abstract type AbstractParticleSource

Type of a particle source that creates primary particles for Geant4 simulations.
Currently defined are [`MonoenergeticSource`](@ref) and [`IsotopeSource`](@ref).
"""
abstract type AbstractParticleSource end


"""
	struct MonoenergeticSource <: AbstractParticleSource

Particle source which emits monoenergetic particles,
in either isotropic, directed ray or directed cone emission.

## Fields
* `particle_type::String`: Particle type emitted by the source.
* `energy::RealQuantity`: Energy of the particles emitted by the source.
* `position::CartesianPoint`: Location of the source
* `direction::CartesianVector`: Direction in which the source emits the particles.
* `opening_angle::AngleQuantity`: Opening angle in case of directed cone emission.

## Examples
```julia 
# Isotropic emission when no direction is passed
m1 = MonoenergeticSource("gamma", 2.615u"MeV", CartesianPoint(0.04, 0, 0.05))

# Directed ray emission when no opening_angle is passed
m2 = MonoenergeticSource("gamma", 2.615u"MeV", CartesianPoint(0.04, 0, 0.05), CartesianVector(-1, 0, 0))

# Directed cone emission when opening_angle is passed
m3 = MonoenergeticSource("gamma", 2.615u"MeV", CartesianPoint(0.04, 0, 0.05), CartesianVector(-1, 0, 0), 10u"°")
```

See also [`IsotopeSource`](@ref).
"""
struct MonoenergeticSource <: AbstractParticleSource
	particle_type::String
	energy::RealQuantity
	position::CartesianPoint
	direction::CartesianVector
	opening_angle::AngleQuantity
	function MonoenergeticSource(particle_type::String, energy::RealQuantity, position::CartesianPoint, direction::CartesianVector, opening_angle::AngleQuantity = 0u"°")
		iszero(direction) && throw(ArgumentError("Direction of the source cannot be zero."))
		(opening_angle < 0u"°" || opening_angle > 180u"°") && throw(ArgumentError("Opening angle of source must be between 0 and 180°."))
		new(particle_type, energy, position, direction, opening_angle)
	end
end

# MonoenergeticSource with isotropic emission <=> opening_angle == 180u"°"
function MonoenergeticSource(particle_type::String, energy::RealQuantity, position::CartesianPoint)
	MonoenergeticSource(particle_type, energy, position, CartesianVector(0,0,1), 180u"°")
end


"""
	struct IsotopeSource <: AbstractParticleSource

Particle source which emits particles based on a defined decay chain,
in either isotropic, directed ray or directed cone emission.

## Fields
* `Z::Int32`: Atomic number of the mother nucleus.
* `A::Int32`: Mass number of the mother nucleus.
* `ionCharge::Float64`: Charge of the mother nucleus.
* `excitEnergy::Float64`: Excitation energy of the mother nucleus.
* `position::CartesianPoint`: Location of the source
* `direction::CartesianVector`: Direction in which the source emits the particles.
* `opening_angle::AngleQuantity`: Opening angle in case of directed cone emission.

## Examples
```julia
# Isotropic emission when no direction is passed
i1 = IsotopeSource(82, 212, 0, 0, CartesianPoint(0.04, 0, 0.05))
   
# Directed ray emission when no opening_angle is passed
i2 = IsotopeSource(82, 212, 0, 0, CartesianPoint(0.04, 0, 0.05), CartesianVector(-1, 0, 0))

# Directed cone emission when opening_angle is passed
i3 = IsotopeSource(82, 212, 0, 0, CartesianPoint(0.04, 0, 0.05), CartesianVector(-1, 0, 0), 10u"°")
```

See also [`MonoenergeticSource`](@ref).
"""
struct IsotopeSource <: AbstractParticleSource
	Z::Int32
	A::Int32
	ionCharge::Float64
	excitEnergy::Float64
	position::CartesianPoint
	direction::CartesianVector
	opening_angle::AngleQuantity
	function IsotopeSource(Z::Integer, A::Integer, ionCharge::Real, excitEnergy::Real, position::CartesianPoint, direction::CartesianVector, opening_angle::AngleQuantity = 0u"°")
		iszero(direction) && throw(ArgumentError("Direction of the source cannot be zero."))
		(opening_angle < 0u"°" || opening_angle > 180u"°") && throw(ArgumentError("Opening angle of source must be between 0 and 180°."))
		new(Int32(Z), Int32(A), Float64(ionCharge), Float64(excitEnergy), position, direction, opening_angle)
	end
end

# IsotopeSource with isotropic emission <=> opening_angle == 180u"°"
function IsotopeSource(Z::Integer, A::Integer, ionCharge::Real, excitEnergy::Real, position::CartesianPoint)
	IsotopeSource(Z, A, ionCharge, excitEnergy, position, CartesianVector(0,0,1), 180u"°")
end

export AbstractParticleSource, MonoenergeticSource, IsotopeSource