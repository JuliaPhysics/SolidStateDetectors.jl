# This file is a part of SolidStateDetectors.jl, licensed under the MIT License (MIT).





function Geometry(T::DataType, geometries::Vector, inputunit::Unitful.Units)#::AbstractGeometry{T} where {T <: AbstractFloat}
    return [Geometry(T, Val{Symbol(geometry["type"])}(), geometry, inputunit) for geometry in geometries]
end