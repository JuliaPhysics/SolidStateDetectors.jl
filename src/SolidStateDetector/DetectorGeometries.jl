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

function get_charge_density(sc::Semiconductor{T}, pt::AbstractCoordinatePoint{T})::T where {T <: SSDFloat}
    get_impurity_density(sc.impurity_density_model, pt) * elementary_charge
end
function get_charge_density(p::Passive{T}, pt::AbstractCoordinatePoint{T})::T where {T <: SSDFloat}
    get_charge_density(p.charge_density_model, pt)
end

function get_ρ_and_ϵ(pt::AbstractCoordinatePoint{T}, det::SolidStateDetector{T}, medium::NamedTuple = material_properties[materials["vacuum"]])::Tuple{T, T, T} where {T <: SSDFloat}
    ρ_semiconductor::T = 0
    q_eff_fix::T = 0
    ϵ::T = medium.ϵ_r
    if pt in det.semiconductor
        ρ_semiconductor = get_charge_density(det.semiconductor, pt) 
        ϵ = det.semiconductor.material.ϵ_r
    elseif !ismissing(det.passives) && in(pt, det.passives)
        for passive in det.passives
            if pt in passive
                q_eff_fix = get_charge_density(passive, pt)
                ϵ = passive.material.ϵ_r
                break
            end
        end
    end
    return ρ_semiconductor, ϵ, q_eff_fix
end