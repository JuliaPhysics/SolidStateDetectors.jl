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

function Tube{T}(dict::Dict{Any, Any}, inputunit::Unitful.Units)::Tube{T} where {T <: SSDFloat}
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


function get_important_points(t::Tube{T}, ::Val{:r})::Vector{T} where {T <: SSDFloat}
    return T[t.r_interval.left, t.r_interval.right]
end

function get_important_points(t::Tube{T}, ::Val{:φ})::Vector{T} where {T <: SSDFloat}
    return T[t.φ_interval.left, t.φ_interval.right]
end

function get_important_points(t::Tube{T}, ::Val{:z})::Vector{T} where {T <: SSDFloat}
    return T[t.z_interval.left, t.z_interval.right]
end

function get_important_points(t::Tube{T}, ::Val{:x})::Vector{T} where {T <: SSDFloat}
    @warn "Not yet implemented"
    return T[]
end

function get_important_points(t::Tube{T}, ::Val{:y})::Vector{T} where {T <: SSDFloat}
    @warn "Not yet implemented"
    return T[]
end

function sample(c::Tube{T}, stepsize::Vector{T}) where T
    samples = CylindricalPoint[]
    for r in c.r_interval.left:stepsize[1]: c.r_interval.right
        for φ in c.φ_interval.left:stepsize[2]: c.φ_interval.right
            for z in c.z_interval.left:stepsize[3]: c.z_interval.right
                push!(samples, CylindricalPoint{T}(r,φ,z))
            end
        end
    end
    return samples
end
