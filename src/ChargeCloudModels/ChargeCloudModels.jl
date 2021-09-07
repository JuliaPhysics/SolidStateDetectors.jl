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
"""
struct NBodyChargeCloud{T <: SSDFloat, N, SH} <: AbstractChargeCloud
    points::Vector{CartesianPoint{T}}
    energies::Vector{T} # in units of eV
    shell_structure::SH
end


"""
    NBodyChargeCloud(center::CartesianPoint{T}, energy::T, particle_type::Type{PT} = Gamma; kwargs...)

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
NBodyChargeCloud(center, energy, number_of_shells = 3, shell_structure = SolidStateDetectors.Icosahedron)
```

!!! note
    Using values with units for `energy` requires the package [Unitful.jl](https://github.com/PainterQubits/Unitful.jl).
"""
function NBodyChargeCloud(center::CartesianPoint{T}, energy::RealQuantity, particle_type::Type{PT} = Gamma;
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

function create_regular_sphere(center::CartesianPoint{T}, N::Integer, R::T)::Vector{CartesianPoint{T}} where {T <: SSDFloat}
    
    # Source: https://www.cmu.edu/biolphys/deserno/pdf/sphere_equi.pdf
    
    # There might be fluctuations of the order 10 points or 1% (whatever is larger)
    points = Vector{CartesianPoint{T}}(undef, N + max(10, round(Int,1.01*N)))
    Ncount = 0
    a = 4π / N
    Mθ = round(π/sqrt(a))
    dθ = π/Mθ
    for m in 0:Mθ-1
        θ = π * (m + 0.5)/Mθ
        Mφ = round(2π * sin(θ)*dθ/a)
        for n in 0:Mφ-1
            φ = 2π*n/Mφ
            Ncount += 1
            points[Ncount] = CartesianPoint{T}(R*cos(φ)*sin(θ), R*sin(φ)*sin(θ), R*cos(θ))
        end
    end
    
    return points[1:Ncount]
end

function NBodyChargeCloud(center::CartesianPoint{T}, energy::RealQuantity, N::Integer, particle_type::Type{PT} = Gamma;
        radius::T = radius_guess(T(to_internal_units(energy)), particle_type), number_of_shells::Int = 2,
    )::NBodyChargeCloud{T, number_of_shells, NTuple{number_of_shells, Int}} where {T, PT <: ParticleType}
    
    points::Vector{CartesianPoint{T}} = CartesianPoint{T}[center]
    energies::Vector{T} = T[1]
    n = Int[]
    
    n_shell = 1
    expected_points_in_shells::Int = div(4 * (4^(number_of_shells) - 1), 3)
    while n_shell <= number_of_shells
        shell = create_regular_sphere(center, round(Int, N * (4^n_shell)/expected_points_in_shells), radius * n_shell)
        points = vcat(points, shell)
        push!(n, length(shell))
        energies = vcat(energies, [exp(-n_shell^2 / 2) for i in 1:length(shell)])
        n_shell += 1
    end
    
    return NBodyChargeCloud{T, number_of_shells, NTuple{number_of_shells, Int}}( points, energies./sum(energies) * T(to_internal_units(energy)), Tuple(n))
end