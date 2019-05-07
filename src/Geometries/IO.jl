# This file is a part of SolidStateDetectors.jl, licensed under the MIT License (MIT).





function Geometries(T::DataType, geometries::Vector, inputunit::Unitful.Units)::Vector{AbstractGeometry}
    return [Geometry(T, Val{Symbol(geometry["type"])}(), Dict{Any,Any}(geometry)::Union{Dict{String,Any},Dict{Any,Any}}, inputunit) for geometry in geometries]
end
