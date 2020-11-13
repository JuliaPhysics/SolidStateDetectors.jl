struct Tube{T} <: AbstractVolumePrimitive{T, 3} ## Only upright Tubes at the moment
    r_interval::AbstractInterval{T}
    φ_interval::AbstractInterval{T}
    z_interval::AbstractInterval{T}
    translate::Union{CartesianVector{T},Missing}
end


function Tube{T}(dict::Dict{Any, Any}, inputunit_dict::Dict{String,Unitful.Units})::Tube{T} where {T <: SSDFloat}

    z_offset::T = 0.0
    if haskey(dict, "translate")
        translate = CartesianVector{T}( haskey(dict["translate"],"x") ? geom_round(ustrip(uconvert(u"m", T(dict["translate"]["x"]) * inputunit_dict["length"] ))) : T(0),
                                        haskey(dict["translate"],"y") ? geom_round(ustrip(uconvert(u"m", T(dict["translate"]["y"]) * inputunit_dict["length"] ))) : T(0),
                                        haskey(dict["translate"],"z") ? geom_round(ustrip(uconvert(u"m", T(dict["translate"]["z"]) * inputunit_dict["length"] ))) : T(0))
        if translate[1] == T(0.0) && translate[2] == T(0.0)
            translate = missing
            z_offset = geom_round(ustrip(uconvert(u"m", T(dict["translate"]["z"]) * inputunit_dict["length"] )))
        end
    else
        translate = missing
    end

    if haskey(dict,"h")
         zStart, zStop = z_offset, geom_round(z_offset + ustrip(uconvert(u"m", T(dict["h"]) * inputunit_dict["length"] )))
    elseif haskey(dict,"z")
         zStart, zStop =  geom_round(z_offset + ustrip(uconvert(u"m", T(dict["z"]["from"]) * inputunit_dict["length"] ))),  geom_round( z_offset + ustrip(uconvert(u"m", T(dict["z"]["to"]) * inputunit_dict["length"])))
    else
        @warn "please specify a height of the Tube 'h'"
    end
    
    φ_interval =  if haskey(dict, "phi") 
        Interval(geom_round(T(ustrip(uconvert(u"rad", T(dict["phi"]["from"]) * inputunit_dict["angle"])))), 
                 geom_round(T(ustrip(uconvert(u"rad", T(dict["phi"]["to"]) * inputunit_dict["angle"])))))                      
    else
        Interval(T(0), geom_round(T(2π)))
    end

    r_interval = if haskey(dict["r"], "from")
        Interval(geom_round(ustrip(uconvert(u"m", T(dict["r"]["from"]) * inputunit_dict["length"] ))), geom_round(ustrip(uconvert(u"m", T(dict["r"]["to"]) * inputunit_dict["length"]))))
    else
        Interval(T(0), geom_round(ustrip(uconvert(u"m", T(dict["r"]) * inputunit_dict["length"]))))
    end

    return Tube{T}(
        r_interval,
        φ_interval,
        Interval(zStart, zStop),
        translate
    )
end

function Geometry(T::DataType, t::Val{:tube}, dict::Dict{Union{Any,String}, Any},inputunit_dict::Dict{String,Unitful.Units})
    return Tube{T}(Dict{Any,Any}(dict), inputunit_dict)
end

function in(point::CartesianPoint{T}, tube::Tube{T})::Bool where {T <: SSDFloat}
    (ismissing(tube.translate) || tube.translate == CartesianVector{T}(0.0,0.0,0.0)) ? nothing : point -= tube.translate
    point = convert(CylindricalPoint,point)
    return (point.r in tube.r_interval && point.φ in tube.φ_interval && point.z in tube.z_interval)
end

function in(point::CylindricalPoint{T}, tube::Tube{T})::Bool where {T <: SSDFloat}
    (ismissing(tube.translate) || tube.translate == CartesianVector{T}(0.0,0.0,0.0)) ? nothing  : point = CylindricalPoint(CartesianPoint(point)-tube.translate)
    return (point.r in tube.r_interval && point.φ in tube.φ_interval && point.z in tube.z_interval)
end


function get_important_points(t::Tube{T}, ::Val{:r})::Vector{T} where {T <: SSDFloat}
    return geom_round.(T[t.r_interval.left, t.r_interval.right])
end

function get_important_points(t::Tube{T}, ::Val{:φ})::Vector{T} where {T <: SSDFloat}
    return geom_round.(T[t.φ_interval.left, t.φ_interval.right])
end

function get_important_points(t::Tube{T}, ::Val{:z})::Vector{T} where {T <: SSDFloat}
    return geom_round.(T[t.z_interval.left, t.z_interval.right])
end

function get_important_points(t::Tube{T}, ::Val{:x})::Vector{T} where {T <: SSDFloat}
    return geom_round.(T[-t.r_interval.left, -t.r_interval.right, t.r_interval.left, t.r_interval.right])
end

function get_important_points(t::Tube{T}, ::Val{:y})::Vector{T} where {T <: SSDFloat}
    return geom_round.(T[-t.r_interval.right, -t.r_interval.left, t.r_interval.left, t.r_interval.right])
end

function sample(c::Tube{T}, stepsize::Vector{T}) where {T <: SSDFloat}
    samples = CylindricalPoint[]
    for r in c.r_interval.left:stepsize[1]: c.r_interval.right
        for φ in c.φ_interval.left:stepsize[2]: c.φ_interval.right
            for z in c.z_interval.left:stepsize[3]: c.z_interval.right
                push!(samples, CylindricalPoint{T}(r,φ,z))
            end
        end
    end
    ismissing(c.translate) ? nothing : samples = map(x -> CylindricalPoint(CartesianPoint(x) + c.translate), samples)
    return samples
end

function (+)(t::Tube{T}, translate::Union{CartesianVector{T},Missing})::Tube{T} where {T <: SSDFloat}
    if ismissing(translate)
        return t
    elseif ismissing(t.translate)
        return Tube(t.r_interval, t.φ_interval, t.z_interval, translate)
    else
        return Tube(t.r_interval, t.φ_interval, t.z_interval, t.translate + translate)
    end
 end
