# """
#     abstract type AbstractConfig{T <: SSDFloat} end
# 
# Supertype of all detector/world/object configs.
# 
# User defined geometries must be subtype of `AbstractConfig{T}`.
# 
# There are a few functions which must be defined for a user config, e.g. `struct UserConfig{T} <: AbstractConfig{T}`:
# 
# For cylindrical grids:
# 
# * in(pt::CylindricalPoint{T}, config::UserConfig{T})::Bool where {T <: SSDFloat}
# * Grid(config::UserConfig{T})::CylindricalGrid{T} where {T <: SSDFloat}
# * get\\_ρ\\_and\\_ϵ(pt::CylindricalPoint{T}, config::UserConfig{T})::Tuple{T, T} where {T <: SSDFloat} 
# * set\\_point\\_types\\_and\\_fixed\\_potentials!(point_types::Array{PointType, 3}, potential::Array{T, 3}, 
#         grid::CylindricalGrid{T}, config::UserConfig{T}; weighting\\_potential\\_contact\\_id::Union{Missing, Int} = missing)::Nothing where {T <: SSDFloat}
# 
# For cartesian grids:
# 
# * in(pt::CartesianPoint{3, T}, config::UserConfig{T})::Bool 
# * Grid(config::UserConfig{T})::CartesianGrid3D{T} where {T <: SSDFloat}
# * get\\_ρ\\_and\\_ϵ(pt::CartesianPoint{3, T}, config::UserConfig{T})::Tuple{T, T} where {T <: SSDFloat} 
# * set\\_point\\_types\\_and\\_fixed\\_potentials!(point_types::Array{PointType, 3}, potential::Array{T, 3}, 
#         grid::CartesianGrid3D{T}, config::UserConfig{T}; weighting\\_potential\\_contact\\_id::Union{Missing, Int} = missing)::Nothing where {T <: SSDFloat}
# 
# """
abstract type AbstractConfig{T <: SSDFloat} end

