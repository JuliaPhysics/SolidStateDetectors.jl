struct Sphere{T} <: AbstractVolumePrimitive{T, 3}
    org::CartesianPoint{T}
    r::T
end

function in(pt::CartesianPoint{T}, s::Sphere{T})::Bool where {T}
    v::CartesianVector{T} = s.org - pt
    return norm(v) <= s.r
end
function in(pt::CylindricalPoint{T}, s::Sphere{T})::Bool where {T}
    return in(CartesianPoint(pt), s)
end

function Sphere{T}(dict::Union{Dict{Any, Any}, Dict{String, Any}}, inputunit_dict::Dict{String,Unitful.Units})::Sphere{T} where {T <: SSDFloat}
    translate_vector = _get_translate_vector(T, dict, inputunit_dict)
    org::CartesianPoint{T} = ismissing(translate_vector) ? CartesianPoint{T}(0, 0, 0) : translate_vector
    r::T = ustrip(uconvert(u"m", dict["r"] * inputunit_dict["length"]))  
    return Sphere{T}(org, r)
end

function Geometry(T::DataType, t::Val{:sphere}, dict::Dict{Any, Any}, inputunit_dict::Dict{String,Unitful.Units})
    return Sphere{T}(dict, inputunit_dict)
end

function (+)(s::Sphere{T}, translate::Union{CartesianVector{T}, Missing})::Sphere{T} where {T <: SSDFloat}
    return ismissing(translate) ? s : Sphere{T}( s.org +  translate, s.r )
end

@recipe function f(s::Sphere{T};) where {T <: SSDFloat}
    @series begin
        pts = []
        label --> ""
        for φ in range(T(0), length = 36, stop = T(2π))
            pt = s.org + CartesianVector{T}( s.r * cos(φ), s.r * sin(φ), 0  )
            push!(pts, pt)
        end
        lines = LineSegment{T, 3, :cartesian}[]
        for i in 1:length(pts)-1
            push!(lines, LineSegment(pts[i+1], pts[i]))
        end
        lines
    end
    @series begin
        label --> ""
        pts = []
        for φ in range(T(0), length = 36, stop = T(2π))
            pt = s.org + CartesianVector{T}( s.r * cos(φ), 0, s.r * sin(φ) )
            push!(pts, pt)
        end
        lines = LineSegment{T, 3, :cartesian}[]
        for i in 1:length(pts)-1
            push!(lines, LineSegment(pts[i+1], pts[i]))
        end
        lines
    end
    @series begin
        pts = []
        label --> ""      
        for φ in range(T(0), length = 36, stop = T(2π))
            pt = s.org + CartesianVector{T}( 0, s.r * cos(φ), s.r * sin(φ) )
            push!(pts, pt)
        end
        lines = LineSegment{T, 3, :cartesian}[]
        for i in 1:length(pts)-1
            push!(lines, LineSegment(pts[i+1], pts[i]))
        end
        lines
    end
end

# For proper grid creation we also need the function get_important_points:
function get_important_points(s::Sphere{T}, ::Val{:r})::Vector{T} where {T <: SSDFloat}
    v::Vector{T} = geom_round.(T[ s.org.x - s.r, s.org.x, s.org.x + s.r, 
                                  s.org.y - s.r, s.org.y, s.org.y + s.r ]) 
    return findall(r -> r >= 0 , v)
end
function get_important_points(s::Sphere{T}, ::Val{:φ})::Vector{T} where {T <: SSDFloat}
    return T[ ]
end
function get_important_points(s::Sphere{T}, ::Val{:z})::Vector{T} where {T <: SSDFloat}
    return geom_round.(T[ s.org.z - s.r, s.org.z, s.org.z + s.r])
end
function get_important_points(s::Sphere{T}, ::Val{:x})::Vector{T} where {T <: SSDFloat}
    return geom_round.(T[ s.org.x - s.r, s.org.x, s.org.x + s.r])
end
function get_important_points(s::Sphere{T}, ::Val{:y})::Vector{T} where {T <: SSDFloat}
    return geom_round.(T[ s.org.y - s.r, s.org.y, s.org.y + s.r])
end

# and a sample function to paint the primitive on the grid (necessary if the object is small)
function sample(s::Sphere{T}, stepsize::Vector{T})::Vector{CartesianPoint{T}}  where {T <: SSDFloat}
    samples::Vector{CartesianPoint{T}} = CartesianPoint{T}[]
    xarr::Vector{T} = get_important_points(s, Val{:x}())
    yarr::Vector{T} = get_important_points(s, Val{:y}())
    zarr::Vector{T} = get_important_points(s, Val{:z}())
    for x in xarr
        for y in yarr
            for z in zarr
                pt::CartesianPoint{T} = CartesianPoint{T}(x, y, z)
                push!(samples, pt)
            end
        end
    end
    return samples
end