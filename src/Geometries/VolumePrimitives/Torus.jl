struct Torus{T} <: AbstractVolumePrimitive{T, 3} ## Only upright Toruss at the moment
    c::T
    a_interval::AbstractInterval{T}
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

    if haskey(dict,"c")
        c = geom_round(ustrip(uconvert(u"m", T(dict["c"]) * inputunit_dict["length"] )))
    else
        @warn "please specify the radius from the center of the hole to the center of the torus tube 'c'"
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

    if haskey(dict,"a")
        if haskey(dict["a"], "from")
            a_interval =  Interval(geom_round(ustrip(uconvert(u"m", T(dict["a"]["from"]) * inputunit_dict["length"] ))), geom_round(ustrip(uconvert(u"m", T(dict["a"]["to"]) * inputunit_dict["length"]))))
        else
            a_interval = Interval(T(0), geom_round(ustrip(uconvert(u"m", T(dict["a"]) * inputunit_dict["length"]))))
        end
    else
        @warn "please specify the radius of torus tube 'a'"
    end

    return Torus{T}(
        c,
        a_interval,
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
    a = sqrt((point.r - Torus.c)^2 + point.z^2)
    θ = acos((point.r - Torus.c)/a)
    if point.z < 0
        θ = 2π - θ
    end
    return (a in Torus.a_interval && point.φ in Torus.φ_interval && θ in Torus.θ_interval)
end

function in(point::CylindricalPoint{T}, Torus::Torus{T})::Bool where {T <: SSDFloat}
    (ismissing(Torus.translate) || Torus.translate == CartesianVector{T}(0.0,0.0,0.0)) ? nothing  : point = CylindricalPoint(CartesianPoint(point)-Torus.translate)
    a = sqrt((point.r - Torus.c)^2 + point.z^2)
    θ = acos((point.r - Torus.c)/a)
    if point.z < 0
        θ = 2π - θ
    end
    return (a in Torus.a_interval && point.φ in Torus.φ_interval && θ in Torus.θ_interval)
end


function get_important_points(t::Torus{T}, ::Val{:a})::Vector{T} where {T <: SSDFloat}
    return geom_round.(T[t.a_interval.left, t.a_interval.right])
end

function get_important_points(t::Torus{T}, ::Val{:r})::Vector{T} where {T <: SSDFloat}
    return geom_round.(T[t.c-t.a_interval.right, t.a_interval.right+t.c])
end

function get_important_points(t::Torus{T}, ::Val{:φ})::Vector{T} where {T <: SSDFloat}
    return geom_round.(T[t.φ_interval.left, t.φ_interval.right])
end

function get_important_points(t::Torus{T}, ::Val{:θ})::Vector{T} where {T <: SSDFloat}
    return geom_round.(T[t.θ_interval.left, t.θ_interval.right])
end

function get_important_points(t::Torus{T}, ::Val{:x})::Vector{T} where {T <: SSDFloat}
    return geom_round.(T[-(t.a_interval.right+t.c), -(t.c-t.a_interval.right), (t.c-t.a_interval.right), (t.a_interval.right+t.c)])
end

function get_important_points(t::Torus{T}, ::Val{:y})::Vector{T} where {T <: SSDFloat}
    return geom_round.(T[-(t.a_interval.right+t.c), -(t.c-t.a_interval.right), (t.c-t.a_interval.right), (t.a_interval.right+t.c)])
end

function get_important_points(t::Torus{T}, ::Val{:z})::Vector{T} where {T <: SSDFloat}
    return geom_round.(T[-t.a_interval.right, t.a_interval.right])
end

function sample(t::Torus{T}, stepsize::Vector{T}) where {T <: SSDFloat}
    samples = CylindricalPoint[]
    for a in t.a_interval.left:stepsize[1]: t.a_interval.right
        for φ in t.φ_interval.left:stepsize[2]: t.φ_interval.right
            for θ in t.θ_interval.left:stepsize[3]: t.θ_interval.right
                push!(samples, CylindricalPoint{T}(t.c+a*cos(θ),φ,a*sin(θ)))
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
        return Torus(t.c, t.a_interval, t.φ_interval, t.θ_interval, translate)
    else
        return Torus(t.c, t.a_interval, t.φ_interval, t.θ_interval, t.translate + translate)
    end
 end
