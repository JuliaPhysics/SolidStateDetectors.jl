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

function point_type(det::SolidStateDetector{T}, grid::Grid{T, 3}, p::CylindricalPoint{T})::Tuple{UInt8, Int, CartesianVector{T}} where {T <: SSDFloat}
    surface_normal::CartesianVector{T} = CartesianVector{T}(0, 0, 0) # need undef version for this
    for contact in det.contacts
        if in(searchsortednearest(grid, p), contact) || in(searchsortednearest(grid, p), contact) || p in contact || p in contact 
            return CD_ELECTRODE::UInt8, contact.id, surface_normal
        end
    end
    on_surface, surface_normal = is_surface_point_and_normal_vector(det, p)
    if on_surface
        for contact in det.contacts
            if in(searchsortednearest(grid, p), contact) #&& abs(sum(sp[2])) > 1
                return CD_ELECTRODE::UInt8, contact.id, surface_normal
            else
                return CD_FLOATING_BOUNDARY::UInt8, -1, surface_normal
            end
        end
    elseif !(p in det)
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
function point_type(det::SolidStateDetector{T}, grid::Grid{T, 3}, p::CartesianPoint{T})::Tuple{UInt8, Int, CartesianVector{T}} where {T <: SSDFloat}
    surface_normal::CartesianVector{T} = CartesianVector{T}(0, 0, 0) # need undef version for this
    for contact in det.contacts
        if p in contact
            return CD_ELECTRODE::UInt8, contact.id, surface_normal
        end
    end
    on_surface, surface_normal = is_surface_point_and_normal_vector(det, p) # surface_normal::CartesianVector{T}
    if on_surface
        for contact in det.contacts
            if in(searchsortednearest(grid, p), contact)
                return CD_ELECTRODE::UInt8, contact.id, surface_normal
            else
                return CD_FLOATING_BOUNDARY::UInt8, -1, surface_normal
            end
        end
    elseif !(p in det)
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

function is_surface_point_and_normal_vector(det::SolidStateDetector{T}, p::CylindricalPoint{T})::Tuple{Bool, CartesianVector{T}} where {T <: SSDFloat}
    if !(p in det) # contacts are already checked in 
        return false, CartesianPoint{T}(0, 0, 0)
    end
    n::MVector{3,T} = @MVector T[0, 0, 0]
    look_around::Vector{Bool} = [   CylindricalPoint{T}(prevfloat(p.r), p.φ, p.z) in det,
                                    CylindricalPoint{T}(nextfloat(p.r), p.φ, p.z) in det,
                                    CylindricalPoint{T}(p.r, prevfloat(p.φ), p.z) in det,
                                    CylindricalPoint{T}(p.r, nextfloat(p.φ), p.z) in det,
                                    CylindricalPoint{T}(p.r, p.φ, prevfloat(p.z)) in det,
                                    CylindricalPoint{T}(p.r, p.φ, nextfloat(p.z)) in det]
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
        Rα::SMatrix{3,3,T} = @SArray([cos(p.φ) -1*sin(p.φ) 0;sin(p.φ) cos(p.φ) 0;0 0 1])
        # return true, geom_round(CartesianVector{T}((Rα * n)...))
        return true, CartesianVector{T}((Rα * n)...)
    end
end
function is_surface_point_and_normal_vector(det::SolidStateDetector{T}, p::CartesianPoint{T})::Tuple{Bool, CartesianVector{T}} where T
    if !(p in det) 
        return false, CartesianPoint{T}(0, 0, 0)
    end
    n::MVector{3,T} = @MVector T[0.0,0.0,0.0]
    look_around::Vector{Bool} = [   CartesianPoint{T}(prevfloat(p.x), p.y, p.z) in det,
                                    CartesianPoint{T}(nextfloat(p.x), p.y, p.z) in det,
                                    CartesianPoint{T}(p.x, prevfloat(p.y), p.z) in det,
                                    CartesianPoint{T}(p.x, nextfloat(p.y), p.z) in det,
                                    CartesianPoint{T}(p.x, p.y, prevfloat(p.z)) in det,
                                    CartesianPoint{T}(p.x, p.y, nextfloat(p.z)) in det]
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
        for ep in det.passives
            if pt in ep
                q_eff_fix = get_charge_density(ep, pt)
                ϵ = ep.material.ϵ_r
                break
            end
        end
    end
    return ρ_semiconductor, ϵ, q_eff_fix
end