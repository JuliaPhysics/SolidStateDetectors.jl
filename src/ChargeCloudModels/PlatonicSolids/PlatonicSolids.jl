struct PlatonicSolid{N, T} <: AbstractChargeCloud
    points::SVector{N, CartesianPoint{T}}
end

function get_shell_structure_points(structure::NTuple{A,Tuple{Int,T,T}}, center::CartesianPoint{T}, length::T = T(1))::SVector where {A,T <: SSDFloat}
    N::Int = sum(map(p -> p[1], structure))
    points::Vector{CartesianPoint{T}} = Vector{CartesianPoint{T}}(undef, N)
    i::Int = 1
    for s in structure
        n::Int = s[1]
        z::T = length * s[2]
        ρ::T = length * sqrt(1 - s[2]^2)
        for φ in range(s[3], step = T(2π/n), length = n)
            sinφ::T, cosφ::T = sincos(φ)
            points[i] = center + CartesianVector{T}(ρ*cosφ,ρ*sinφ,z)
            i+=1
        end
    end
    SVector{N, CartesianPoint{T}}(points)
end

get_vertices(::Type{PlatonicSolid{N, T}}) where {N, T} = N


"""
    struct PointCharge{T} <: AbstractChargeCloud

Struct which defines a single point-like charge carrier.

## Fields
* `points::SVector{1, CartesianPoint{T}}`: Position of the charge carrier, saved as single entry of a `Vector`.

## Example
```julia
center = CartesianPoint{T}(0,0,0)
pc = PointCharge(center) # Constructor: creates a PointCharge around (0,0,0)
pc.points                # Array that only contains the center point
```
"""
const PointCharge{T} = PlatonicSolid{1, T}
PointCharge(center::Vector{CartesianPoint{T}}, args...) where {T} = PointCharge{T}(center)
PointCharge(center::CartesianPoint{T}, args...) where {T} = PointCharge{T}([center])


const Tetrahedron{T} = PlatonicSolid{4, T}
function Tetrahedron(center::CartesianPoint{T}, length::Real = T(1))::Tetrahedron{T} where {T}
    structure::NTuple{2, Tuple{Int,T,T}} = ((1,T(1),T(0)), (3,T(-1/3),T(0)))
    PlatonicSolid(get_shell_structure_points(structure, center, T(length)))
end

const Hexahedron{T} = PlatonicSolid{8, T}
function Hexahedron(center::CartesianPoint{T}, length::Real = T(1))::Hexahedron{T} where {T}
    structure::NTuple{4, Tuple{Int,T,T}} = ((1,T(1),T(0)),(3,T(1/3),T(0)),(3,T(-1/3),T(π/3)),(1,T(-1),T(0)))
    PlatonicSolid(get_shell_structure_points(structure, center, T(length)))
end

const Octahedron{T} = PlatonicSolid{6, T}
function Octahedron(center::CartesianPoint{T}, length::Real = T(1))::Octahedron{T} where {T}
    structure::NTuple{3, Tuple{Int,T,T}} = ((1,T(1),T(0)), (4,T(0),T(0)), (1,T(-1),T(0)))
    PlatonicSolid(get_shell_structure_points(structure, center, T(length)))
end

const Icosahedron{T} = PlatonicSolid{12, T}
function Icosahedron(center::CartesianPoint{T}, length::Real = T(1))::Icosahedron{T} where {T}
    structure::NTuple{4, Tuple{Int,T,T}} = ((1,T(1),T(0)), (5,T(sqrt(1/5)),T(0)), (5,T(-sqrt(1/5)),T(π/5)), (1,T(-1),T(0)))
    PlatonicSolid(get_shell_structure_points(structure, center, T(length)))
end

const Dodecahedron{T} = PlatonicSolid{20, T}
function Dodecahedron(center::CartesianPoint{T}, length::Real = T(1))::Dodecahedron{T} where {T}
    zmax::T = sqrt((5 + 2*sqrt(5))/15)
    tmp::T = 1 - 2*(1 - zmax^2) * sin(π/5)^2
    zmin::T = zmax * tmp - sqrt(1 - zmax^2) * sqrt(1 - tmp^2)
    structure::NTuple{4, Tuple{Int,T,T}} = ((5,zmax,T(0)),(5,zmin,T(0)),(5,-zmin,T(π/5)),(5,-zmax,T(π/5)))
    PlatonicSolid(get_shell_structure_points(structure, center, T(length)))
end
