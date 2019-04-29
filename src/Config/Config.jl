"""
    abstract type AbstractConfig{T <: SSDFloat} end

Supertype of all detector/world/object configs.

User defined geometries must be subtype of `AbstractConfig{T}`.

There are a few functions which must be defined for a user config, e.g. `struct UserConfig{T} <: AbstractConfig{T}`:

For cylindrical grids:

* in(pt::cylindrical{T}, config::UserConfig{T})::Bool where {T <: SSDFloat}
* Grid(config::UserConfig{T})::Grid{T, 3, :cylindrical} where {T <: SSDFloat}
* get\\_ρ\\_and\\_ϵ(pt::cylindrical{T}, config::UserConfig{T})::Tuple{T, T} where {T <: SSDFloat} 
* set\\_pointtypes\\_and\\_fixed\\_potentials!(pointtypes::Array{PointType, 3}, potential::Array{T, 3}, 
        grid::Grid{T, 3, :cylindrical}, config::UserConfig{T}; weighting\\_potential\\_channel\\_idx::Union{Missing, Int} = missing)::Nothing where {T <: SSDFloat}

For cartesian grids:

* in(pt::StaticVector{3, T}, config::UserConfig{T})::Bool 
* Grid(config::UserConfig{T})::Grid{T, 3, :cartesian} where {T <: SSDFloat}
* get\\_ρ\\_and\\_ϵ(pt::StaticVector{3, T}, config::UserConfig{T})::Tuple{T, T} where {T <: SSDFloat} 
* set\\_pointtypes\\_and\\_fixed\\_potentials!(pointtypes::Array{PointType, 3}, potential::Array{T, 3}, 
        grid::Grid{T, 3, :cartesian}, config::UserConfig{T}; weighting\\_potential\\_channel\\_idx::Union{Missing, Int} = missing)::Nothing where {T <: SSDFloat}

"""
abstract type AbstractConfig{T <: SSDFloat} end

