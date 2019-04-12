
struct ConeMantle{T} <: AbstractGeometry{T, 3} ## Only upright Cones at the moment, Convention: counterclockwise \alpha \beta γ; γ is the 90 deg angle, 
    cone::Cone{T}
end

function Geometry(T::DataType, t::Val{:ConeMantle}, dict::Dict{Any, Any}, inputunit::Unitful.Units)
    return ConeMantle{T}(Cone{T}(dict, inputunit))
end

function in(point::CartesianPoint{T}, cm::ConeMantle)::Bool where T
    convert(CylindricalPoint{T}, point)
    point in cm 
end

function in(point::CylindricalPoint{T}, cm::ConeMantle{T})::Bool where T
    r = get_diagonal_r_from_z(cm.cone, point.z)
    if point.r == r && point.φ in cm.cone.φ_interval && point.z in cm.cone.z_interval
        return true
    else
        return false
    end
end

function in(point::CylindricalPoint{T}, cm::ConeMantle{T}, rs::Vector{T})::Bool where T
    analytical_r = get_diagonal_r_from_z(cm.cone,point.z)
    r = rs[searchsortednearest(rs,analytical_r)]
    if point.r == r && point.φ in cm.cone.φ_interval && point.z in cm.cone.z_interval
        return true
    else
        return false
    end
end



function get_important_points(c::ConeMantle{T}, s)::Vector{T} where {T <: SSDFloat}
    get_important_points(c.cone, s)
end

