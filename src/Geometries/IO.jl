# This file is a part of SolidStateDetectors.jl, licensed under the MIT License (MIT).

function Geometry(T::DataType, geometries::Vector, inputunit::Unitful.Units)#::AbstractGeometry{T} where {T <: SSDFloat}
    return [Geometry(T, Val{Symbol(geometry["type"])}(), Dict{Any,Any}(geometry)::Union{Dict{String,Any},Dict{Any,Any}}, inputunit) for geometry in geometries]
end
