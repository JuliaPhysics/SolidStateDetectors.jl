orientations = Dict("tl" => :top_left,
                    "lt" => :top_left,
                    "top_left" => :top_left,
                    "left_top" => :top_left,
                    "bl" =>:bottom_left,
                    "lb" =>:bottom_left,
                    "bottom_left" => :bottom_left,
                    "left_bottom" => :bottom_left,
                    "tr" => :top_right,
                    "rt" => :top_right,
                    "top_right" => :top_right,
                    "right_top" => :top_right,
                    "br" => :bottom_right,
                    "rb" => :bottom_right,
                    "bottom_right" => :bottom_right,
                    "right_bottom" => :bottom_right)




struct Cone{T} <: AbstractGeometry{T, 3} ## Only upright Cones at the moment, Convention: counterclockwise \alpha \beta γ; γ is the 90 deg angle, 
    org::AbstractCoordinatePoint{T}
    hierarchy::Int
    polarity::Symbol
    orientation::Symbol
    r_interval::AbstractInterval{T}
    φ_interval::AbstractInterval{T}
    z_interval::AbstractInterval{T}
    α::T  #Convention: counterclockwise α β γ; γ is the 90 deg angle
    β::T
end



function Cone{T}(org::AbstractCoordinatePoint{T}, hierarchy::Int, polarity::String, orientation::String, rStart::T, rStop::T, φStart::T, φStop::T, zStart::T, zStop::T) where {T}
    o = orientations[orientation]
    o == :top_right || o == :bottom_left ? α = atan(abs(zStop-zStart)/abs(rStop-rStart)) : α = atan((rStop-rStart)/abs(zStop-zStart))
    β = π/2 -α
    pol = pols[polarity]
    if pol == :pos
        Cone{T}(org, hierarchy, pol, o, ClosedInterval{T}(rStart, rStop) ,ClosedInterval{T}(deg2rad(φStart) ,deg2rad(φStop)),ClosedInterval{T}(zStart, zStop), α, β)
    elseif pol == :neg
        if deg2rad(φStop)%2π == deg2rad(φStart)
            Cone{T}(org, hierarchy, pol, o, OpenInterval{T}(rStart, rStop) ,OpenInterval{T}(deg2rad(φStart) ,deg2rad(φStop)),OpenInterval{T}(zStart, zStop), α, β)
        else
            Cone{T}(org, hierarchy, pol, o, OpenInterval{T}(rStart, rStop) ,ClosedInterval{T}(deg2rad(φStart) ,deg2rad(φStop)),OpenInterval{T}(zStart, zStop), α, β)
        end
    end
end

function Cone(hierarchy::Int, polarity::String, orientation::String ,rStart::T, rStop::T, φStart::T, φStop::T, zStart::T, zStop::T) where T
    Cone{T}(CartesianPoint{T}(0.0,0.0,0.0), hierarchy, polarity, orientation, rStart, rStop, φStart, φStop, zStart, zStop)
end

function Cone{T}(hierarchy::Int, polarity::String, orientation::String, r_interval::AbstractInterval{T}, φ_interval::AbstractInterval{T}, z_interval::AbstractInterval{T}) where T
    Cone(hierarchy, polarity, orientation, r_interval.left, r_interval.right, φ_interval.left, φ_interval.right, z_interval.left,z_interval.right)
end

function in(point::CylindricalPoint{T}, cone::Cone{T}) where T
    if cone.orientation == :top_right || cone.orientation == :bottom_right
        cone.polarity == :neg ? x = (point.r > get_diagonal_r_from_z(cone, point.z) && point.r < cone.r_interval.right) : x = (point.r >= get_diagonal_r_from_z(cone, point.z) && point.r <= cone.r_interval.right)
        if x && point.φ in cone.φ_interval && point.z in cone.z_interval
            return true 
        end
    elseif cone.orientation == :top_left || cone.orientation == :bottom_left
        cone.polarity == :neg ? x = (point.r >  cone.r_interval.left && point.r < get_diagonal_r_from_z(cone, point.z)) : x = (point.r >=  cone.r_interval.left && point.r <= get_diagonal_r_from_z(cone, point.z))
        if x  && point.φ in cone.φ_interval && point.z in cone.z_interval
            return true 
        end
    end
    return false
end

function in(point::CartesianPoint{T}, cone::Cone{T}) where T
    point = convert(CylindricalPoint,point)
    point in cone
end

function get_diagonal_r_from_z(cone::Cone{T}, z::T) where T
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

function Cone{T}(dict::Dict{Any, Any}, inputunit::Unitful.Units)::Cone{T} where {T <: SSDFloat}
    haskey(dict, "hierarchy") ? h = dict["hierarchy"] : h = 1
    haskey(dict, "pol") ? polarity = dict["pol"] : polarity = "positive"
    return Cone{T}(h,
        polarity,
        dict["orientation"],
        Interval(geom_round(ustrip(uconvert(u"m", T(dict["rStart"]) * inputunit ))), geom_round(ustrip(uconvert(u"m", T(dict["rStop"]) * inputunit)))),
        Interval(geom_round(T(dict["phiStart"])), geom_round(T(dict["phiStop"]))),
        Interval(geom_round(ustrip(uconvert(u"m", T(dict["zStart"]) * inputunit ))), geom_round(ustrip(uconvert(u"m", T(dict["zStop"]) * inputunit))))  )
end

function Geometry(T::DataType, t::Val{:Cone}, dict::Dict{Any, Any}, inputunit::Unitful.Units)
    return Cone{T}(dict, inputunit)
end


function get_important_points(c::Cone{T}, ::Val{:r})::Vector{T} where {T <: SSDFloat}
    return T[c.r_interval.left, c.r_interval.right] 
end

function get_important_points(c::Cone{T}, ::Val{:φ})::Vector{T} where {T <: SSDFloat}
    return T[c.φ_interval.left, c.φ_interval.right]
end

function get_important_points(c::Cone{T}, ::Val{:z})::Vector{T} where {T <: SSDFloat}
    return T[c.z_interval.left, c.z_interval.right]
end

function get_important_points(c::Cone{T}, ::Val{:x})::Vector{T} where {T <: SSDFloat}
    @warn "Not yet implemented"
    return T[]
end
function get_important_points(c::Cone{T}, ::Val{:y})::Vector{T} where {T <: SSDFloat}
    @warn "Not yet implemented"
    return T[]
end

