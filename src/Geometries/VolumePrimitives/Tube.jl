struct Tube{T} <: AbstractGeometry{T, 3} ## Only upright Tubes at the moment
    org::CartesianPoint{T}
    hierarchy::Int
    polarity::Symbol
    r_interval::AbstractInterval{T}
    φ_interval::AbstractInterval{T}
    z_interval::AbstractInterval{T}
end

function Tube{T}(org::CartesianPoint{T}, hierarchy::Int, polarity::String, rStart::T,rStop::T,φStart::T,φStop::T,zStart::T,zStop::T) where {T}
    pol = pols[polarity]
    if pol == :neg
        if deg2rad(φStop)%T(2π) == deg2rad(φStart)
            Tube{T}(org, hierarchy, pol, OpenInterval{T}(rStart, rStop), ClosedInterval{T}(deg2rad(φStart) ,deg2rad(φStop)),OpenInterval{T}(zStart ,zStop))
        else
            Tube{T}(org, hierarchy, pol, OpenInterval{T}(rStart, rStop), OpenInterval{T}(deg2rad(φStart) ,deg2rad(φStop)),OpenInterval{T}(zStart ,zStop))
        end
    else
        Tube{T}(org, hierarchy, pol, ClosedInterval{T}(rStart, rStop), ClosedInterval{T}(deg2rad(φStart) ,deg2rad(φStop)),ClosedInterval{T}(zStart ,zStop))
    end
end

function Tube{T}(hierarchy::Int, polarity::String, rStart::T,rStop::T,φStart::T,φStop::T,zStart::T,zStop::T) where {T}
    Tube{T}(CartesianPoint{T}(0.0,0.0,0.0), hierarchy, polarity, rStart,rStop,φStart,φStop,zStart,zStop)
end

function Tube{T}(
    hierarchy::Int,
    polarity::String,
    r_interval::AbstractInterval{T},
    φ_interval::AbstractInterval{T},
    z_interval::AbstractInterval{T}) where T
    return Tube{T}(hierarchy, polarity, r_interval.left, r_interval.right, φ_interval.left, φ_interval.right, z_interval.left, z_interval.right)
end

function Tube{T}(dict::Dict{Any, Any}, inputunit::Unitful.Units)::Tube{T} where {T <: AbstractFloat}
    haskey(dict, "hierarchy") ? h = dict["hierarchy"] : h = 1
    haskey(dict, "pol") ? polarity = dict["pol"] : polarity = "positive"
    return Tube{T}(h,
    polarity,
    Interval(geom_round(ustrip(uconvert(u"m", T(dict["rStart"]) * inputunit ))), geom_round(ustrip(uconvert(u"m", T(dict["rStop"]) * inputunit)))),
    Interval(geom_round(T(dict["phiStart"])), geom_round(T(dict["phiStop"]))),
    Interval(geom_round(ustrip(uconvert(u"m", T(dict["zStart"]) * inputunit ))), geom_round(ustrip(uconvert(u"m", T(dict["zStop"]) * inputunit))))  )
end

function Geometry(T::DataType, t::Val{:Tube}, dict::Dict{Union{Any,String}, Any}, inputunit::Unitful.Units)
    return Tube{T}(Dict{Any,Any}(dict), inputunit)
end

function in(point::CartesianPoint{T}, tube::Tube{T}) where T
    point = convert(CylindricalPoint,point)
    if point.r in tube.r_interval && point.φ in tube.φ_interval && point.z in tube.z_interval
        return true
    end
    return false
end

function in(point::CylindricalPoint{T}, tube::Tube{T}) where T
    if point.r in tube.r_interval && point.φ in tube.φ_interval && point.z in tube.z_interval
        return true
    end
    return false
end


function get_important_points(t::Tube{T})::NTuple{3, Vector{T}} where {T <: AbstractFloat}
    v1::Vector{T} = T[t.r_interval.left, t.r_interval.right] #[g.x[1], g.x[2]]
    v2::Vector{T} = T[t.φ_interval.left, t.φ_interval.right]
    v3::Vector{T} = T[t.z_interval.left, t.z_interval.right]
    return v1, v2, v3
end