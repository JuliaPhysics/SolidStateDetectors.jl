include("Object.jl")
include("Passive.jl")
include("Contacts.jl")
include("Semiconductor.jl")
include("VirtualVolumes.jl")
include("SolidStateDetector.jl")

function sample(det::SolidStateDetector{T}, ::Type{Cartesian}, sampling...)::Vector{CartesianPoint{T}} where {T <: SSDFloat}
    imp::Vector{CartesianPoint{T}} = vcat(
        CartesianPoint.(sample(det.semiconductor.geometry, sampling...)),
        [CartesianPoint.(sample(g.geometry, sampling...)) for object in skipmissing((det.contacts, det.passives)) for g in object]...)
    unique!(imp)
end

function sample(det::SolidStateDetector{T}, ::Type{Cylindrical}, sampling...)::Vector{CylindricalPoint{T}} where {T <: SSDFloat}
    imp::Vector{CylindricalPoint{T}} = vcat(
    CylindricalPoint.(sample(det.semiconductor.geometry, sampling...)),
    [CylindricalPoint.(sample(g.geometry, sampling...)) for object in skipmissing((det.contacts, det.passives)) for g in object]...)
    unique!(imp)
end


# go_to_nearest_gridpoint
function searchsortednearest(grid::CylindricalGrid{T}, pt::CylindricalPoint{T})::CylindricalPoint{T} where {T <: SSDFloat}
    idx1::Int = searchsortednearest(grid.axes[1].ticks, pt.r)
    idx2::Int = searchsortednearest(grid.axes[2].ticks, pt.φ)
    idx3::Int = searchsortednearest(grid.axes[3].ticks, pt.z)
    CylindricalPoint{T}(grid.axes[1].ticks[idx1], grid.axes[2].ticks[idx2], grid.axes[3].ticks[idx3])
end
function searchsortednearest(grid::CartesianGrid3D{T}, pt::CartesianPoint{T})::CartesianPoint{T} where {T <: SSDFloat}
    idx1::Int = searchsortednearest(grid.axes[1].ticks, pt.x)
    idx2::Int = searchsortednearest(grid.axes[2].ticks, pt.y)
    idx3::Int = searchsortednearest(grid.axes[3].ticks, pt.z)
    CartesianPoint{T}(grid.axes[1].ticks[idx1], grid.axes[2].ticks[idx2], grid.axes[3].ticks[idx3])
end

get_charge_density(sc::Semiconductor{T}, pt::AbstractCoordinatePoint{T}) where {T <: SSDFloat} = 
    get_impurity_density(sc.impurity_density_model, pt) * elementary_charge

get_charge_density(p::Passive{T}, pt::AbstractCoordinatePoint{T}) where {T <: SSDFloat} =
    get_charge_density(p.charge_density_model, pt)

@inline function get_ρimp_ϵ_ρfix(pt::AbstractCoordinatePoint{T}, obj::Semiconductor)::Tuple{T, T, T} where {T <: SSDFloat}
    return get_charge_density(obj, pt), obj.material.ϵ_r, zero(T)
end
@inline function get_ρimp_ϵ_ρfix(pt::AbstractCoordinatePoint{T}, obj::Passive)::Tuple{T, T, T} where {T <: SSDFloat}
    return zero(T), obj.material.ϵ_r, get_charge_density(obj, pt) 
end