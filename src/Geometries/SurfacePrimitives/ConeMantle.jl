
struct ConeMantle{T} <: AbstractGeometry{T, 3} ## Only upright Cones at the moment, Convention: counterclockwise \alpha \beta γ; γ is the 90 deg angle,
    hierarchy::Int
    cone::Cone{T}
end

function Geometry(T::DataType, t::Val{:ConeMantle}, dict::Dict{Any, Any}, inputunit::Unitful.Units)
    return ConeMantle{T}( Cone{T}(dict, inputunit).hierarchy, Cone{T}(dict, inputunit) )
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


function get_r(cm::ConeMantle)
    return cm.cone.r_interval.left, cm.cone.r_interval.right
end

function get_φ(cm::ConeMantle)
    return cm.cone.φ_interval.left, cm.cone.φ_interval.right
end

function get_z(cm::ConeMantle)
    return cm.cone.z_interval.left, cm.cone.z_interval.right
end
