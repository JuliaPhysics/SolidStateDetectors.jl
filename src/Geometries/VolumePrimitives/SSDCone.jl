struct SSDCone{T} <: AbstractGeometry{T, 3} ## Only upright Cones at the moment, Convention: counterclockwise \alpha \beta γ; γ is the 90 deg angle,
    hierarchy::Int
    rStart1::T
    rStop1::T
    rStart2::T
    rStop2::T
    φStart::T
    φStop::T
    zStart::T
    zStop::T

    translate::SVector{3,T}
    rotate::SMatrix{3,3,T}
    # α1::T
    # α2::T
    # α::T  #Convention: counterclockwise α β γ; γ is the 90 deg angle
    # β::T
end



function SSDCone{T}(hierarchy::Int, rStart1::T, rStop1::T, rStart2::T, rStop2::T, φStart::T, φStop::T, zStart::T, zStop::T, translate::SVector{3,T}, rotX::T, rotY::T, rotZ::T) where {T}
    rotationMatrix::SMatrix{3,3,T} = SMatrix{3,3,T}([1 0 0;
                                                    0 cos(rotX) -sin(rotX);
                                                    0 sin(rotX) cos(rotX)]) *
                                    SMatrix{3,3,T}([cos(rotY) 0 sin(rotY);
                                                    0 1 0;
                                                    -sin(rotY) 0 cos(rotY)]) *
                                     SMatrix{3,3,T}([cos(rotZ) -sin(rotZ) 0;
                                                    sin(rotZ) cos(rotZ) 0;
                                                    0 0 1])
    # a⃗1::SVector{2,T} = [rStop1 - rStart1, 0.0]
    # b⃗1::SVector{2,T} =  [rStart2, zStop] - [rStart1,zStart]
    # a⃗2::SVector{2,T} = -a⃗1
    # b⃗2::SVector{2,T} =  [rStop2, zStop] - [rStop1,zStart]
    # α1::T = acos( dot(a⃗1,b⃗1) / ( norm(a⃗1) * norm(b⃗1) ) )
    # α2::T = acos( dot(a⃗2,b⃗2) / ( norm(a⃗2) * norm(b⃗2) ) )

    SSDCone{T}(hierarchy, rStart1, rStop1, rStart2, rStop2, φStart, φStop, zStart, zStop, translate, rotationMatrix)
end

function in(point::CylindricalPoint{T}, cone::SSDCone{T}) where T
    if ( point.z in ClosedInterval{T}(cone.zStart,cone.zStop) ) && ( point.φ in ClosedInterval{T}(cone.φStart,cone.φStop) ) && ( point.r in ClosedInterval{T}(minimum([cone.rStart1,cone.rStart2,cone.rStop1,cone.rStop2]), maximum([cone.rStart1,cone.rStart2,cone.rStop1,cone.rStop2])) )

        r1,r2 = get_intersection_rs_for_given_z(point.z,cone)
        if point.r in ClosedInterval{T}(r1,r2)
            return true
        else
            return false
        end

    else
        return false
    end
end

function get_intersection_rs_for_given_z(z::T, c::SSDCone{T}) where {T}
    r1::T = ( (c.rStart1*c.zStop - c.zStart*c.rStart2) * (T(0.0) - T(1.0)) - (c.rStart1 - c.rStart2) * (T(0.0)*z - z*T(1.0)) ) / ( (c.rStart1 - c.rStart2) * (z - z) - (c.zStart - c.zStop) * (T(0.0) - T(1.0)) )
    r2::T = ( (c.rStop1*c.zStop - c.zStart*c.rStop2) * (T(0.0) - T(1.0)) - (c.rStop1 - c.rStop2) * (T(0.0)*z - z*T(1.0)) ) / ( (c.rStop1 - c.rStop2) * (z - z) - (c.zStart - c.zStop) * (T(0.0) - T(1.0)) )
    r1,r2
end

function in(point::CartesianPoint{T}, cone::SSDCone{T}) where T
    point = convert(CylindricalPoint,point)
    point in cone
end

function get_diagonal_r_from_z(cone::SSDCone{T}, z::T) where T
    if cone.orientation == :top_right
        r = tan(cone.β)*(z-cone.z_interval.left)
        return geom_round(T(cone.r_interval.right - r))
    elseif cone.orientation == :top_left
        r = tan(cone.α)*(z-cone.z_interval.left)
        return geom_round(T(cone.r_interval.left + r))
    elseif cone.orientation == :bottom_right
        r = tan(cone.α)*(z-cone.z_interval.left)
        return geom_round(T(cone.r_interval.left + r))
    elseif cone.orientation == :bottom_left
        r =  tan(cone.β)*(cone.z_interval.right -z)
        return geom_round(T(cone.r_interval.left + r))
    end
end

function SSDCone{T}(dict::Union{Dict{Any, Any},Dict{String,Any}}, inputunit::Unitful.Units)::SSDCone{T} where {T <: AbstractFloat}
    haskey(dict, "hierarchy") ? h::Int = dict["hierarchy"] : h = 1
    haskey(dict, "translate") ? translate::SVector{3,T} = dict["translate"] : translate = SVector{3,T}(0,0,0)
    haskey(dict, "rotX") ? rotX::T = deg2rad(dict["rotX"]) : rotX = T(0.0)
    haskey(dict, "rotY") ? rotY::T = deg2rad(dict["rotY"]) : rotY = T(0.0)
    haskey(dict, "rotZ") ? rotZ::T = deg2rad(dict["rotZ"]) : rotZ = T(0.0)
    return SSDCone{T}(h,
        geom_round(ustrip(uconvert(u"m", T(dict["rStart1"]) * inputunit ))),
        geom_round(ustrip(uconvert(u"m", T(dict["rStop1"]) * inputunit))),
        geom_round(ustrip(uconvert(u"m", T(dict["rStart2"]) * inputunit ))),
        geom_round(ustrip(uconvert(u"m", T(dict["rStop2"]) * inputunit))),
        geom_round(T(deg2rad(dict["phiStart"]))),
        geom_round(T(deg2rad(dict["phiStop"]))),
        geom_round(ustrip(uconvert(u"m", T(dict["zStart"]) * inputunit ))),
        geom_round(ustrip(uconvert(u"m", T(dict["zStop"]) * inputunit))),
        translate, rotX, rotY, rotZ)
end

function Geometry(T::DataType, t::Val{:SSDCone}, dict::Dict{Any, Any}, inputunit::Unitful.Units)
    return SSDCone{T}(dict, inputunit)
end


function get_important_points(c::SSDCone{T})::NTuple{3, Vector{T}} where {T <: AbstractFloat}
    v1::Vector{T} = T[c.r_interval.left, c.r_interval.right] #[g.x[1], g.x[2]]
    v2::Vector{T} = T[c.φ_interval.left, c.φ_interval.right]
    v3::Vector{T} = T[c.z_interval.left, c.z_interval.right]
    return v1, v2, v3
end