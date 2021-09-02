abstract type AbstractChargeCloud end

include("PlatonicSolids/PlatonicSolids.jl")
include("ParticleTypes.jl")

"""
    struct NBodyChargeCloud{T, N, SH} <: AbstractChargeCloud

Struct which defines a charge cloud consisting of multiple point-like charge carriers,
initially distributed around a given center.

## Parametric Types 
* `T`: Precision type.
* `N`: Number of shells.
* `SH`: Geometry of the shells (should be a `PlatonicSolid`).

## Fields
* `points::Vector{CartesianPoint{T}}`: Positions of the charge carriers that are part of the charge cloud.
* `energies::Vector{T}`: Energies of the respective charge carriers, in the same order as `points`.
* `shell_structure::SH`: Initial geometry of the charge carriers around the `center` point, relevant for plotting.

See also [`create_charge_cloud`](@ref).
"""
struct NBodyChargeCloud{T <: SSDFloat, N, SH} <: AbstractChargeCloud
    points::Vector{CartesianPoint{T}}
    energies::Vector{T} # in units of eV
    shell_structure::SH
end


"""
    create_charge_cloud(center::CartesianPoint{T}, energy::T, particle_type::Type{PT} = Gamma; kwargs...)

Returns an [`NBodyChargeCloud`](@ref) for a given `energy` deposition at a position that defines the `center` of the charge cloud.

## Arguments
* `center::CartesianPoint{T}`: Center position of the [`NBodyChargeCloud`](@ref).
* `energy::RealQuantity`: Deposited energy with units. If no units are given, the value is parsed in units of eV.
* `particle_type`: [`ParticleType`](@ref) of the particle that deposited the energy. Default is `Gamma`.

## Keywords
* `radius::T`: Estimate for the radius of the [`NBodyChargeCloud`](@ref). Default is determined from `particle_type` via `radius_guess`.
* `number_of_shells::Int`: Number of shells around the `center` point. Default is `2`.
* `shell_structure`: Geometry with which the charges are distributed in the shells. Default is `Dodecahedron`.

## Example
```julia
center = CartesianPoint{T}([0,0,0])
energy = 1460u"keV"
create_charge_cloud(center, energy, number_of_shells = 3, shell_structure = SolidStateDetectors.Icosahedron)
```

!!! note
    Using values with units for `energy` requires the package [Unitful.jl](https://github.com/PainterQubits/Unitful.jl).

See also [`NBodyChargeCloud`](@ref).
"""
function create_charge_cloud(center::CartesianPoint{T}, energy::RealQuantity, particle_type::Type{PT} = Gamma;
        radius::T = radius_guess(T(to_internal_units(energy)), particle_type), number_of_shells::Int = 2, shell_structure = Dodecahedron
    )::NBodyChargeCloud{T, number_of_shells, typeof(shell_structure{T})} where {T, PT <: ParticleType}
    
    points::Vector{CartesianPoint{T}} = CartesianPoint{T}[center]
    energies::Vector{T} = T[1]
    
    n_shell = 1
    while n_shell <= number_of_shells
        shell = shell_structure(center, n_shell * radius).points
        points = vcat(points, shell)
        energies = vcat(energies, [exp(-n_shell^2 / 2) for i in 1:length(shell)])
        n_shell += 1
    end
    
    return NBodyChargeCloud{T, number_of_shells, typeof(shell_structure{T})}( points, energies./sum(energies) * T(to_internal_units(energy)), shell_structure{T})
end