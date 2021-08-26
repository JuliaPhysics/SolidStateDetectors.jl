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

function point_type(det::SolidStateDetector{T}, grid::Grid{T, 3}, pt::CylindricalPoint{T})::Tuple{UInt8, Int, CartesianVector{T}} where {T <: SSDFloat}
    surface_normal::CartesianVector{T} = CartesianVector{T}(0, 0, 0) # need undef version for this
    for contact in det.contacts
        if in(searchsortednearest(grid, pt), contact) || in(searchsortednearest(grid, pt), contact) || pt in contact || pt in contact 
            return CD_ELECTRODE::UInt8, contact.id, surface_normal
        end
    end
    on_surface, surface_normal = is_surface_point_and_normal_vector(det, pt)
    if on_surface
        for contact in det.contacts
            if in(searchsortednearest(grid, pt), contact) #&& abs(sum(sp[2])) > 1
                return CD_ELECTRODE::UInt8, contact.id, surface_normal
            else
                return CD_FLOATING_BOUNDARY::UInt8, -1, surface_normal
            end
        end
    elseif !(pt in det)
        return CD_OUTSIDE::UInt8, -1, surface_normal
    else
        return CD_BULK::UInt8, -1, surface_normal
    end
end


# Point types for charge drift
const CD_ELECTRODE = 0x00
const CD_OUTSIDE = 0x01
const CD_BULK = 0x02
const CD_FLOATING_BOUNDARY = 0x04 # not 0x03, so that one could use bit operations here...

# """
#     For charge drift...
# """
function point_type(det::SolidStateDetector{T}, grid::Grid{T, 3}, pt::CartesianPoint{T})::Tuple{UInt8, Int, CartesianVector{T}} where {T <: SSDFloat}
    surface_normal::CartesianVector{T} = CartesianVector{T}(0, 0, 0) # need undef version for this
    for contact in det.contacts
        if pt in contact
            return CD_ELECTRODE::UInt8, contact.id, surface_normal
        end
    end
    on_surface, surface_normal = is_surface_point_and_normal_vector(det, pt) # surface_normal::CartesianVector{T}
    if on_surface
        for contact in det.contacts
            if in(searchsortednearest(grid, pt), contact)
                return CD_ELECTRODE::UInt8, contact.id, surface_normal
            else
                return CD_FLOATING_BOUNDARY::UInt8, -1, surface_normal
            end
        end
    elseif !(pt in det)
        return CD_OUTSIDE::UInt8, -1, surface_normal
    else
        return CD_BULK::UInt8, -1, surface_normal
    end
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

function is_surface_point_and_normal_vector(det::SolidStateDetector{T}, pt::CylindricalPoint{T})::Tuple{Bool, CartesianVector{T}} where {T <: SSDFloat}
    if !(pt in det) # contacts are already checked in 
        return false, CartesianPoint{T}(0, 0, 0)
    end
    n::MVector{3,T} = @MVector T[0, 0, 0]
    look_around::Vector{Bool} = [   CylindricalPoint{T}(prevfloat(pt.r), pt.φ, pt.z) in det,
                                    CylindricalPoint{T}(nextfloat(pt.r), pt.φ, pt.z) in det,
                                    CylindricalPoint{T}(pt.r, prevfloat(pt.φ), pt.z) in det,
                                    CylindricalPoint{T}(pt.r, nextfloat(pt.φ), pt.z) in det,
                                    CylindricalPoint{T}(pt.r, pt.φ, prevfloat(pt.z)) in det,
                                    CylindricalPoint{T}(pt.r, pt.φ, nextfloat(pt.z)) in det]
    if all(look_around)
        return false , CartesianPoint{T}(n...)
    else
        if !look_around[1] n[1] -= 1 end
        if !look_around[2] n[1] += 1 end
        if !look_around[3] n[2] -= 1 end
        if !look_around[4] n[2] += 1 end
        if !look_around[5] n[3] -= 1 end
        if !look_around[6] n[3] += 1 end
        # println(look_around , " " , n)
        Rα::SMatrix{3,3,T} = @SArray([cos(pt.φ) -1*sin(pt.φ) 0;sin(pt.φ) cos(pt.φ) 0;0 0 1])
        # return true, geom_round(CartesianVector{T}((Rα * n)...))
        return true, CartesianVector{T}((Rα * n)...)
    end
end
function is_surface_point_and_normal_vector(det::SolidStateDetector{T}, pt::CartesianPoint{T})::Tuple{Bool, CartesianVector{T}} where T
    if !(pt in det) 
        return false, CartesianPoint{T}(0, 0, 0)
    end
    n::MVector{3,T} = @MVector T[0.0,0.0,0.0]
    look_around::Vector{Bool} = [   CartesianPoint{T}(prevfloat(pt.x), pt.y, pt.z) in det,
                                    CartesianPoint{T}(nextfloat(pt.x), pt.y, pt.z) in det,
                                    CartesianPoint{T}(pt.x, prevfloat(pt.y), pt.z) in det,
                                    CartesianPoint{T}(pt.x, nextfloat(pt.y), pt.z) in det,
                                    CartesianPoint{T}(pt.x, pt.y, prevfloat(pt.z)) in det,
                                    CartesianPoint{T}(pt.x, pt.y, nextfloat(pt.z)) in det]
    if all(look_around)
        return false, n
    else
        if !look_around[1] n[1] -= 1 end
        if !look_around[2] n[1] += 1 end
        if !look_around[3] n[2] -= 1 end
        if !look_around[4] n[2] += 1 end
        if !look_around[5] n[3] -= 1 end
        if !look_around[6] n[3] += 1 end
        return true, n
    end
end


get_charge_density(sc::Semiconductor{T}, pt::AbstractCoordinatePoint{T}) where {T <: SSDFloat} = 
    get_impurity_density(sc.impurity_density_model, pt) * elementary_charge

get_charge_density(p::Passive{T}, pt::AbstractCoordinatePoint{T}) where {T <: SSDFloat} =
    get_charge_density(p.charge_density_model, pt)

@inline function get_ρ_and_ϵ(pt::AbstractCoordinatePoint{T}, obj::Semiconductor)::Tuple{T, T, T} where {T <: SSDFloat}
    q_eff_imp = get_charge_density(obj, pt) 
    ϵ = obj.material.ϵ_r
    return q_eff_imp, ϵ, zero(T)
end
@inline function get_ρ_and_ϵ(pt::AbstractCoordinatePoint{T}, obj::Passive)::Tuple{T, T, T} where {T <: SSDFloat}
    q_eff_fix = get_charge_density(obj, pt) 
    ϵ = obj.material.ϵ_r
    return zero(T), ϵ, q_eff_fix
end