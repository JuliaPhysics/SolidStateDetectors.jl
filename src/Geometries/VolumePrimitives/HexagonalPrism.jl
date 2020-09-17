struct HexagonalPrism{T} <: AbstractVolumePrimitive{T, 3} ## Only upright hexagons at the moment
    a::T #side length/ maximal radius
    h::T #total height
    translate::Union{CartesianPoint{T}, Missing}#origin at middle of hexagonal prism
end

# You also have to implement the function to obtain the primitive from a config file (so an dic)
# You also should provide a example config file containing this new primitive
function HexagonalPrism{T}(dict::Union{Dict{Any, Any}, Dict{String, Any}}, inputunit_dict::Dict{String,Unitful.Units})::HexagonalPrism{T} where {T <: SSDFloat}
    # ... parse values from dict to NewPrimitive{T}(...)
    if haskey(dict, "translate")
        translate::CartesianPoint{T} = CartesianPoint{T}(
            haskey(dict["translate"],"x") ? geom_round(ustrip(uconvert(u"m", T(dict["translate"]["x"]) * inputunit_dict["length"] ))) : T(0),
            haskey(dict["translate"],"y") ? geom_round(ustrip(uconvert(u"m", T(dict["translate"]["y"]) * inputunit_dict["length"] ))) : T(0),
            haskey(dict["translate"],"z") ? geom_round(ustrip(uconvert(u"m", T(dict["translate"]["z"]) * inputunit_dict["length"] ))) : T(0))
    else
        translate = missing
    end
    a::T = ustrip(uconvert(u"m", dict["a"] * inputunit_dict["length"]))
    h::T = ustrip(uconvert(u"m", dict["h"] * inputunit_dict["length"]))
    return HexagonalPrism{T}(a, h, translate)
end

# Auxiliary function for point in a triangle check
function tri_area(p1::CartesianPoint{T}, p2::CartesianPoint{T}, p3::CartesianPoint{T})::T where {T}
    return (p1.x - p3.x) * (p2.y - p3.y) - (p2.x - p3.x) * (p1.y - p3.y)
end

# is a point inside the hexagonal prism?
function in(pt::CartesianPoint{T}, hp::HexagonalPrism{T})::Bool where {T}
    pt = ismissing(hp.translate) ? pt : pt .- hp.translate
    cpt = CylindricalPoint(pt) #shift pt to prism's frame
    minRad = sqrt(3)/2*hp.a
    # Use hexagonal symmetry and do the ol' point in a triangle test, algebraically
    angle = π/6 - max(0, mod(cpt.φ, π/3)-π/6)
    pt_mod = CartesianPoint(CylindricalPoint{T}(cpt.r, angle, cpt.z))

    return abs(pt_mod.z) <= hp.h/2 && # Within the height
            cpt.r <= hp.a && # within outer radius
            if cpt.r > minRad # Inside inner and outer radius
                v1 = CartesianPoint{T}(0, 0, 0)
                v2 = CartesianPoint{T}(minRad, 0, 0)
                v3 = CartesianPoint{T}(minRad, hp.a/2, 0)
                b1 = tri_area(pt_mod, v1, v2) < 0.
                b2 = tri_area(pt_mod, v2, v3) < 0.
                b3 = tri_area(pt_mod, v3, v1) < 0.
                ((b1 == b2) && (b2 == b3))
            else # Inside inner radius
                true
            end
end

function in(pt::CylindricalPoint{T}, hp::HexagonalPrism{T})::Bool where {T}
    return in(CartesianPoint(pt), hp)
end

function Geometry(T::DataType, t::Val{:HexagonalPrism}, dict::Dict{Any, Any}, inputunit_dict::Dict{String,Unitful.Units})
    return HexagonalPrism{T}(dict, inputunit_dict)
end

# add a (+) method to shift the primitive
function (+)(hp::HexagonalPrism{T}, translate::Union{CartesianVector{T}, Missing})::HexagonalPrism{T} where {T <: SSDFloat}
    if ismissing(translate)
        return hp
    elseif ismissing(hp.translate)
        return HexagonalPrism(hp.a, hp.h, translate)
    else
        return HexagonalPrism(hp.a, hp.h, hp.translate + translate)
    end
end

# Also a plot recipe for this new primitive should be provided:
@recipe function f(hp::HexagonalPrism{T}) where {T <: SSDFloat}
    label --> "HexagonalPrism"
    @series begin
        pts_top = []
        pts_bottom = []

        #find all vertices, this loop has been tested and works
        for φ in [0,deg2rad(60), deg2rad(120), deg2rad(180), deg2rad(240), deg2rad(300)]
            pt_top = hp.translate .+ CartesianPoint{T}(hp.a * cos(φ), hp.a * sin(φ), hp.h/2)
            push!(pts_top, pt_top)
            pt_bottom = hp.translate .+ CartesianPoint{T}(hp.a * cos(φ), hp.a * sin(φ), -hp.h/2)
            push!(pts_bottom, pt_bottom)
        end

        #create Linesegments connecting the vertices
        lines = LineSegment{T, 3, :cartesian}[]
        N = length(pts_top)
        for i in 1:N
            push!(lines, LineSegment(pts_top[i%N+1], pts_top[i])) #top hexagon
            push!(lines, LineSegment(pts_bottom[i%N+1], pts_bottom[i])) #bottom hexagon
            push!(lines, LineSegment(pts_bottom[i], pts_top[i])) #lines connecting hexagons (missing one line)
        end
        lines
    end
end

# For proper grid creation we also need the function get_important_points:
function get_important_points(hp::HexagonalPrism{T}, ::Val{:r})::Vector{T} where {T <: SSDFloat}
    return geom_round.(T[sqrt(3)/2*hp.a, hp.a])
end

function get_important_points(hp::HexagonalPrism{T}, ::Val{:φ})::Vector{T} where {T <: SSDFloat}
    return T[]
end
function get_important_points(hp::HexagonalPrism{T}, ::Val{:z})::Vector{T} where {T <: SSDFloat}
    return geom_round.(T[-hp.h/2, hp.h/2])
end
function get_important_points(hp::HexagonalPrism{T}, ::Val{:x})::Vector{T} where {T <: SSDFloat}
    return geom_round.(T[-hp.a, -sqrt(3)/2*hp.a, sqrt(3)/2*hp.a, hp.a])
end
function get_important_points(hp::HexagonalPrism{T}, ::Val{:y})::Vector{T} where {T <: SSDFloat}
    return geom_round.(T[-hp.a, -sqrt(3)/2*hp.a, sqrt(3)/2*hp.a, hp.a])
end

#and a sample function to paint the primitive on the grid (necessary if the object is small)
function sample(hp::HexagonalPrism{T}, stepsize::Vector{T})  where T
    minRad = sqrt(3)/2*hp.a
    samples = CartesianPoint{T}[]
    for x in -hp.a : stepsize[1] : hp.a
        for y in -hp.a : stepsize[2]: hp.a
            for z in -hp.h/2 : stepsize[3] : hp.h/2
                p = CartesianPoint{T}(x, y, z)
                if p.x <= minRad && p.y <= minRad
                    push!(samples, p)
                elseif in(p, hp)
                    push!(samples, p)
                end
            end
        end
    end
    ismissing(hp.translate) ? nothing : samples = map(x -> x + hp.translate, samples)
    return samples
end
