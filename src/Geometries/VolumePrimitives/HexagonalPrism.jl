struct HexagonalPrism{T} <: AbstractVolumePrimitive{T, 3} ## Only upright hexagons at the moment
    rInner::T
    rOuter::T
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
    rInner, rOuter = if haskey(dict["r"], "from")
        geom_round(ustrip(uconvert(u"m", T(dict["r"]["from"]) * inputunit_dict["length"] ))), geom_round(ustrip(uconvert(u"m", T(dict["r"]["to"]) * inputunit_dict["length"])))
    else
        geom_round(0), geom_round(ustrip(uconvert(u"m", T(dict["r"]["to"]) * inputunit_dict["length"])))
    end
    h::T = ustrip(uconvert(u"m", dict["h"] * inputunit_dict["length"]))
    return HexagonalPrism{T}(rInner, rOuter, h, translate)
end

# Auxiliary function for point in a triangle check
function tri_area(p1::CartesianPoint{T}, p2::CartesianPoint{T}, p3::CartesianPoint{T})::T where {T}
    return (p1.x - p3.x) * (p2.y - p3.y) - (p2.x - p3.x) * (p1.y - p3.y)
end

# is a point inside the hexagonal prism?
function in(pt::CartesianPoint{T}, hp::HexagonalPrism{T})::Bool where {T}
    pt = ismissing(hp.translate) ? pt : pt .- hp.translate
    cpt = CylindricalPoint(pt) #shift pt to prism's frame
    minRadInner = sqrt(3)/2*hp.rInner
    minRadOuter = sqrt(3)/2*hp.rOuter
    # Use hexagonal symmetry and do the ol' point in a triangle test, algebraically
    angle = π/6 - max(0, mod(cpt.φ, π/3)-π/6)
    pt_mod = CartesianPoint(CylindricalPoint{T}(cpt.r, angle, cpt.z))

    return abs(pt_mod.z) <= hp.h/2 && # Within the height
            cpt.r <= hp.rOuter && # within outer radius
            if cpt.r >= minRadInner # Inside inner and outer radius
                # v1 = CartesianPoint{T}(0, 0, 0)
                # v2 = CartesianPoint{T}(minRad, 0, 0)
                # v3 = CartesianPoint{T}(minRad, hp.a/2, 0)
                # b1 = tri_area(pt_mod, v1, v2) < 0.
                # b2 = tri_area(pt_mod, v2, v3) < 0.
                # b3 = tri_area(pt_mod, v3, v1) < 0.
                # ((b1 == b2) && (b2 == b3))
                ((pt_mod.x >= minRadInner) && (pt_mod.x <= minRadOuter))
            else # Inside inner radius
                false
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
        return HexagonalPrism(hp.rInner, hp.rOuter, hp.h, translate)
    else
        return HexagonalPrism(hp.rInner, hp.rOuter, hp.h, hp.translate + translate)
    end
end

# Also a plot recipe for this new primitive should be provided:
@recipe function f(hp::HexagonalPrism{T}) where {T <: SSDFloat}
    label --> "HexagonalPrism"
    @series begin
        pts_top_outer = []
        pts_bottom_outer = []
        pts_top_inner = []
        pts_bottom_inner = []

        #find all vertices, this loop has been tested and works
        for φ in [0,deg2rad(60), deg2rad(120), deg2rad(180), deg2rad(240), deg2rad(300)]
            pt_top_outer = hp.translate .+ CartesianPoint{T}(hp.rOuter * cos(φ), hp.rOuter * sin(φ), hp.h/2)
            push!(pts_top_outer, pt_top_outer)
            pt_top_inner = hp.translate .+ CartesianPoint{T}(hp.rInner * cos(φ), hp.rInner * sin(φ), hp.h/2)
            push!(pts_top_inner, pt_top_inner)
            pt_bottom_outer = hp.translate .+ CartesianPoint{T}(hp.rOuter * cos(φ), hp.rOuter * sin(φ), -hp.h/2)
            push!(pts_bottom_outer, pt_bottom_outer)
            pt_bottom_inner = hp.translate .+ CartesianPoint{T}(hp.rInner * cos(φ), hp.rInner * sin(φ), -hp.h/2)
            push!(pts_bottom_inner, pt_bottom_inner)
        end

        #create Linesegments connecting the vertices
        lines = LineSegment{T, 3, :cartesian}[]
        N = length(pts_top_outer)
        for i in 1:N
            push!(lines, LineSegment(pts_top_outer[i%N+1], pts_top_outer[i])) #top outer hexagon
            push!(lines, LineSegment(pts_top_inner[i%N+1], pts_top_inner[i])) #top inner hexagon
            push!(lines, LineSegment(pts_bottom_outer[i%N+1], pts_bottom_outer[i])) #bottom outer hexagon
            push!(lines, LineSegment(pts_bottom_inner[i%N+1], pts_bottom_inner[i])) #bottom inner hexagon
            push!(lines, LineSegment(pts_bottom_outer[i], pts_top_outer[i])) #lines connecting outer hexagons
            push!(lines, LineSegment(pts_bottom_inner[i], pts_top_inner[i])) #lines connecting inner hexagons
            push!(lines, LineSegment(pts_top_outer[i], pts_top_inner[i])) #lines connecting hexagons
            push!(lines, LineSegment(pts_bottom_outer[i], pts_bottom_inner[i])) #lines connecting hexagons
        end
        lines
    end
end

# For proper grid creation we also need the function get_important_points:
# Radial
function get_important_points(hp::HexagonalPrism{T}, ::Val{:r})::Vector{T} where {T <: SSDFloat}
    return geom_round.(T[sqrt(3)/2*hp.rInner, hp.rInner, sqrt(3)/2*hp.rOuter, hp.rOuter])
end

# polar angle
function get_important_points(hp::HexagonalPrism{T}, ::Val{:φ})::Vector{T} where {T <: SSDFloat}
    return T[]
end
# Z
function get_important_points(hp::HexagonalPrism{T}, ::Val{:z})::Vector{T} where {T <: SSDFloat}
    return geom_round.(T[-hp.h/2, hp.h/2])
end
# X
function get_important_points(hp::HexagonalPrism{T}, ::Val{:x})::Vector{T} where {T <: SSDFloat}
    return geom_round.(T[-hp.rOuter, -sqrt(3)/2*hp.rOuter, -hp.rInner, -sqrt(3)/2*hp.rInner, sqrt(3)/2*hp.rInner, hp.rInner, sqrt(3)/2*hp.rOuter, hp.rOuter])
end
# Y
function get_important_points(hp::HexagonalPrism{T}, ::Val{:y})::Vector{T} where {T <: SSDFloat}
    return geom_round.(T[-hp.rOuter, -sqrt(3)/2*hp.rOuter, -hp.rInner, -sqrt(3)/2*hp.rInner, sqrt(3)/2*hp.rInner, hp.rInner, sqrt(3)/2*hp.rOuter, hp.rOuter])
end

#and a sample function to paint the primitive on the grid (necessary if the object is small)
function sample(hp::HexagonalPrism{T}, stepsize::Vector{T})  where T
    minRadInner = sqrt(3)/2*hp.rInner
    minRadOuter = sqrt(3)/2*hp.rOuter
    samples = CartesianPoint{T}[]
    for x in -hp.rOuter : stepsize[1] : hp.rOuter
        for y in -hp.rOuter : stepsize[2]: hp.rOuter
            for z in -hp.h/2 : stepsize[3] : hp.h/2
                p = CartesianPoint{T}(x, y, z)
                if in(samples, p)
                    push!(samples, p)
                end
                # if p.x >= minRad && p.y >= minRad
                #     push!(samples, p)
                # elseif in(p, hp)
                #     push!(samples, p)
                # end
            end
        end
    end
    ismissing(hp.translate) ? nothing : samples = map(x -> x + hp.translate, samples)
    return samples
end
