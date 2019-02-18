"""
    abstract type AbstractConfig{T <: AbstractFloat} end

Supertype of all detector/world/object configs.

User defined geometries must be subtype of `AbstractConfig{T}`.

There are a few functions which must be defined for a user config, e.g. `struct UserConfig{T} <: AbstractConfig{T}`:

For cylindrical grids:

* in(pt::Cylindrical{T}, config::UserConfig{T})::Bool where {T <: AbstractFloat}
* Grid(config::UserConfig{T})::Grid{T, 3, :Cylindrical} where {T <: AbstractFloat}
* get\\_ρ\\_and\\_ϵ(pt::Cylindrical{T}, config::UserConfig{T})::Tuple{T, T} where {T <: AbstractFloat} 
* set\\_pointtypes\\_and\\_fixed\\_potentials!(pointtypes::Array{PointType, 3}, potential::Array{T, 3}, 
        grid::Grid{T, 3, :Cylindrical}, config::UserConfig{T}; weighting\\_potential\\_channel\\_idx::Union{Missing, Int} = missing)::Nothing where {T <: AbstractFloat}

For cartesian grids:

* in(pt::StaticVector{3, T}, config::UserConfig{T})::Bool 
* Grid(config::UserConfig{T})::Grid{T, 3, :Cartesian} where {T <: AbstractFloat}
* get\\_ρ\\_and\\_ϵ(pt::StaticVector{3, T}, config::UserConfig{T})::Tuple{T, T} where {T <: AbstractFloat} 
* set\\_pointtypes\\_and\\_fixed\\_potentials!(pointtypes::Array{PointType, 3}, potential::Array{T, 3}, 
        grid::Grid{T, 3, :Cartesian}, config::UserConfig{T}; weighting\\_potential\\_channel\\_idx::Union{Missing, Int} = missing)::Nothing where {T <: AbstractFloat}

"""
abstract type AbstractConfig{T <: AbstractFloat} end

