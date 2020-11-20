"""
    # struct Torus{T}

r_torus: Radius measured from origin to center of tube in the standard r coordinate.
r_tube_interval: Internal to external radius of cross section of torus tube.
φ_interval: Standard cylindrical coordinate φ
θ_interval: Angle measured from plane of torus, in a perpendicular cross section.
...
"""

struct Torus{T} <: AbstractVolumePrimitive{T, 3} ## Only upright Torus at the moment
    r_torus::T
    r_tube_interval::AbstractInterval{T}
    φ_interval::AbstractInterval{T}
    θ_interval::AbstractInterval{T}
    translate::Union{CartesianVector{T},Missing}
end


function Torus{T}(dict::Dict{Any, Any}, inputunit_dict::Dict{String,Unitful.Units})::Torus{T} where {T <: SSDFloat}

    if haskey(dict, "translate")
        translate = CartesianVector{T}( haskey(dict["translate"],"x") ? geom_round(ustrip(uconvert(u"m", T(dict["translate"]["x"]) * inputunit_dict["length"] ))) : T(0),
                                        haskey(dict["translate"],"y") ? geom_round(ustrip(uconvert(u"m", T(dict["translate"]["y"]) * inputunit_dict["length"] ))) : T(0),
                                        haskey(dict["translate"],"z") ? geom_round(ustrip(uconvert(u"m", T(dict["translate"]["z"]) * inputunit_dict["length"] ))) : T(0))
    else
        translate = missing
    end

    if haskey(dict,"r_torus")
        r_torus = geom_round(ustrip(uconvert(u"m", T(dict["r_torus"]) * inputunit_dict["length"] )))
    else
        @warn "please specify the radius from the center of the hole to the center of the torus tube 'r_torus'"
    end

    θ_interval =  if haskey(dict, "theta")
        Interval(geom_round(T(ustrip(uconvert(u"rad", T(dict["theta"]["from"]) * inputunit_dict["angle"])))),
                 geom_round(T(ustrip(uconvert(u"rad", T(dict["theta"]["to"]) * inputunit_dict["angle"])))))
    else
        Interval(T(0), geom_round(T(2π)))
    end

    φ_interval =  if haskey(dict, "phi")
        Interval(geom_round(T(ustrip(uconvert(u"rad", T(dict["phi"]["from"]) * inputunit_dict["angle"])))),
                 geom_round(T(ustrip(uconvert(u"rad", T(dict["phi"]["to"]) * inputunit_dict["angle"])))))
    else
        Interval(T(0), geom_round(T(2π)))
    end

    if haskey(dict,"r_tube")
        if haskey(dict["r_tube"], "from")
            r_tube_interval =  Interval(geom_round(ustrip(uconvert(u"m", T(dict["r_tube"]["from"]) * inputunit_dict["length"] ))), geom_round(ustrip(uconvert(u"m", T(dict["r_tube"]["to"]) * inputunit_dict["length"]))))
        else
            r_tube_interval = Interval(T(0), geom_round(ustrip(uconvert(u"m", T(dict["r_tube"]) * inputunit_dict["length"]))))
        end
    else
        @warn "please specify the radius of torus tube 'r_tube'"
    end

    return Torus{T}(
        r_torus,
        r_tube_interval,
        φ_interval,
        θ_interval,
        translate
    )
end

function Geometry(T::DataType, t::Val{:torus}, dict::Dict{Union{Any,String}, Any},inputunit_dict::Dict{String,Unitful.Units})
    return Torus{T}(Dict{Any,Any}(dict), inputunit_dict)
end

function in(point::CartesianPoint{T}, Torus::Torus{T})::Bool where {T <: SSDFloat}
    (ismissing(Torus.translate) || Torus.translate == CartesianVector{T}(0.0,0.0,0.0)) ? nothing : point -= Torus.translate
    point = convert(CylindricalPoint,point)
    r_tube = sqrt((point.r - Torus.r_torus)^2 + point.z^2)
    θ = acos((point.r - Torus.r_torus)/r_tube)
    if point.z < 0
        θ = 2π - θ
    end
    return (r_tube in Torus.r_tube_interval && point.φ in Torus.φ_interval && θ in Torus.θ_interval)
end

function in(point::CylindricalPoint{T}, Torus::Torus{T})::Bool where {T <: SSDFloat}
    (ismissing(Torus.translate) || Torus.translate == CartesianVector{T}(0.0,0.0,0.0)) ? nothing  : point = CylindricalPoint(CartesianPoint(point)-Torus.translate)
    r_tube = sqrt((point.r - Torus.r_torus)^2 + point.z^2)
    θ = acos((point.r - Torus.r_torus)/r_tube)
    if point.z < 0
        θ = 2π - θ
    end
    return (r_tube in Torus.r_tube_interval && point.φ in Torus.φ_interval && θ in Torus.θ_interval)
end


function get_important_points(t::Torus{T}, ::Val{:r_tube})::Vector{T} where {T <: SSDFloat}
    return geom_round.(T[t.r_tube_interval.left, t.r_tube_interval.right])
end

function get_important_points(t::Torus{T}, ::Val{:r})::Vector{T} where {T <: SSDFloat}
    return geom_round.(T[t.r_torus-t.r_tube_interval.right, t.r_tube_interval.right+t.r_torus])
end

function get_important_points(t::Torus{T}, ::Val{:φ})::Vector{T} where {T <: SSDFloat}
    return geom_round.(T[t.φ_interval.left, t.φ_interval.right])
end

function get_important_points(t::Torus{T}, ::Val{:θ})::Vector{T} where {T <: SSDFloat}
    return geom_round.(T[t.θ_interval.left, t.θ_interval.right])
end

function get_important_points(t::Torus{T}, ::Val{:x})::Vector{T} where {T <: SSDFloat}
    return geom_round.(T[-(t.r_tube_interval.right+t.r_torus), -(t.r_torus-t.r_tube_interval.right), (t.r_torus-t.r_tube_interval.right), (t.r_tube_interval.right+t.r_torus)])
end

function get_important_points(t::Torus{T}, ::Val{:y})::Vector{T} where {T <: SSDFloat}
    return geom_round.(T[-(t.r_tube_interval.right+t.r_torus), -(t.r_torus-t.r_tube_interval.right), (t.r_torus-t.r_tube_interval.right), (t.r_tube_interval.right+t.r_torus)])
end

function get_important_points(t::Torus{T}, ::Val{:z})::Vector{T} where {T <: SSDFloat}
    return geom_round.(T[-t.r_tube_interval.right, t.r_tube_interval.right])
end

function sample(t::Torus{T}, stepsize::Vector{T}) where {T <: SSDFloat}
    samples = CylindricalPoint[]
    for r_tube in t.r_tube_interval.left:stepsize[1]: t.r_tube_interval.right
        for φ in t.φ_interval.left:stepsize[2]: t.φ_interval.right
            for θ in t.θ_interval.left:stepsize[3]: t.θ_interval.right
                push!(samples, CylindricalPoint{T}(t.r_torus+r_tube*cos(θ),φ,r_tube*sin(θ)))
            end
        end
    end
    ismissing(t.translate) ? nothing : samples = map(x -> CylindricalPoint(CartesianPoint(x) + t.translate), samples)
    return samples
end

function (+)(t::Torus{T}, translate::Union{CartesianVector{T},Missing})::Torus{T} where {T <: SSDFloat}
    if ismissing(translate)
        return t
    elseif ismissing(t.translate)
        return Torus(t.r_torus, t.r_tube_interval, t.φ_interval, t.θ_interval, translate)
    else
        return Torus(t.r_torus, t.r_tube_interval, t.φ_interval, t.θ_interval, t.translate + translate)
    end
 end
