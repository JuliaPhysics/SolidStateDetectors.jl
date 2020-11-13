struct HexagonalPrism{T} <: AbstractVolumePrimitive{T, 3} ## Only upright hexagons at the moment
    rInner::T
    rOuter::T
    h::T #total height
    translate::Union{CartesianVector{T}, Missing} #origin at middle of hexagonal prism
    rotZ::T #rotation around z-axis, counterclockwise when looking from the top
end

# You also have to implement the function to obtain the primitive from a config file (so an dic)
# You also should provide a example config file containing this new primitive
function HexagonalPrism{T}(dict::Union{Dict{Any, Any}, Dict{String, Any}}, inputunit_dict::Dict{String,Unitful.Units})::HexagonalPrism{T} where {T <: SSDFloat}

    if haskey(dict, "translate")
        translate = CartesianPoint{T}(
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
    rotZ::T = haskey(dict, "rotZ") ? geom_round(T(ustrip(uconvert(u"rad", T(dict["rotZ"]) * inputunit_dict["angle"])))) : T(0)
    return HexagonalPrism{T}(rInner, rOuter, h, translate, rotZ)
end

# is a point inside the hexagonal prism?
function in(pt::CartesianPoint{T}, hp::HexagonalPrism{T})::Bool where {T <: SSDFloat}
    pt = ismissing(hp.translate) ? pt : pt - hp.translate
    cpt = CylindricalPoint(pt) #shift pt to prism's frame
    cpt = CylindricalPoint{T}( cpt.r, cpt.φ + hp.rotZ, cpt.z ) #rotate it by rotZ
    minRadInner = sqrt(3)/2*hp.rInner
    minRadOuter = sqrt(3)/2*hp.rOuter
    # Use hexagonal symmetry to sweep the point into a basic triangle
    angle = if (mod(cpt.φ ÷ (π/6), 2) > 0) π/6 - mod(cpt.φ, (π/6)) else mod(cpt.φ, (π/6)) end

    return abs(pt.z) <= hp.h/2 && # Within the height
            cpt.r <= hp.rOuter && # within outer radius
            cpt.r >= minRadInner && # Inside minimal inner radius
            ((cpt.r*cos(angle) >= minRadInner) && (cpt.r*cos(angle) <= minRadOuter))
end

@inline function in(pt::CylindricalPoint{T}, hp::HexagonalPrism{T})::Bool where {T <: SSDFloat}
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
        return HexagonalPrism(hp.rInner, hp.rOuter, hp.h, translate, hp.rotZ)
    else
        return HexagonalPrism(hp.rInner, hp.rOuter, hp.h, hp.translate + translate, hp.rotZ)
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
        for φ in deg2rad.(30:60:330) .- hp.rotZ
            pt_top_outer = CartesianPoint{T}(hp.translate.x + hp.rOuter * cos(φ), hp.translate.y + hp.rOuter * sin(φ), hp.translate.z + hp.h/2)
            push!(pts_top_outer, pt_top_outer)
            pt_top_inner = CartesianPoint{T}(hp.translate.x + hp.rInner * cos(φ), hp.translate.y + hp.rInner * sin(φ), hp.translate.z + hp.h/2)
            push!(pts_top_inner, pt_top_inner)
            pt_bottom_outer = CartesianPoint{T}(hp.translate.x + hp.rOuter * cos(φ), hp.translate.y + hp.rOuter * sin(φ), hp.translate.z -hp.h/2)
            push!(pts_bottom_outer, pt_bottom_outer)
            pt_bottom_inner = CartesianPoint{T}(hp.translate.x + hp.rInner * cos(φ), hp.translate.y + hp.rInner * sin(φ), hp.translate.z -hp.h/2)
            push!(pts_bottom_inner, pt_bottom_inner)
        end

        #create Linesegments connecting the vertices
        lines = LineSegment{T, 3, :cartesian}[]
        N = length(pts_top_outer)
        for i in 1:N
            push!(lines, LineSegment(pts_top_outer[i%N+1], pts_top_outer[i])) #top outer hexagon
            push!(lines, LineSegment(pts_bottom_outer[i%N+1], pts_bottom_outer[i])) #bottom outer hexagon
            push!(lines, LineSegment(pts_bottom_outer[i], pts_top_outer[i])) #lines connecting outer hexagons
            if hp.rInner > 0
                push!(lines, LineSegment(pts_bottom_inner[i%N+1], pts_bottom_inner[i])) #bottom inner hexagon
                push!(lines, LineSegment(pts_top_inner[i%N+1], pts_top_inner[i])) #top inner hexagon
                push!(lines, LineSegment(pts_bottom_inner[i], pts_top_inner[i])) #lines connecting inner hexagons
                push!(lines, LineSegment(pts_top_outer[i], pts_top_inner[i])) #lines connecting hexagons
                push!(lines, LineSegment(pts_bottom_outer[i], pts_bottom_inner[i])) #lines connecting hexagons
            end
        end
        lines
    end
end

# For proper grid creation we also need the function get_important_points:
# Radial
function get_important_points(hp::HexagonalPrism{T}, ::Val{:r})::Vector{T} where {T <: SSDFloat}
    return geom_round.(T[-hp.rOuter, -sqrt(3)/2*hp.rOuter, -hp.rInner, -sqrt(3)/2*hp.rInner, 0., sqrt(3)/2*hp.rInner, hp.rInner, sqrt(3)/2*hp.rOuter, hp.rOuter] .+ sqrt(hp.translate.x^2 + hp.translate.y^2))
end

# polar angle
function get_important_points(hp::HexagonalPrism{T}, ::Val{:φ})::Vector{T} where {T <: SSDFloat}
    return T[]
end
# Z
function get_important_points(hp::HexagonalPrism{T}, ::Val{:z})::Vector{T} where {T <: SSDFloat}
    return geom_round.(T[-hp.h/2, hp.h/2] .+ hp.translate.z)
end
# X
function get_important_points(hp::HexagonalPrism{T}, ::Val{:x})::Vector{T} where {T <: SSDFloat}
    return geom_round.(T[-hp.rOuter, -sqrt(3)/2*hp.rOuter, -hp.rInner, -sqrt(3)/2*hp.rInner, 0., sqrt(3)/2*hp.rInner, hp.rInner, sqrt(3)/2*hp.rOuter, hp.rOuter] .+ hp.translate.x)
end
# Y
function get_important_points(hp::HexagonalPrism{T}, ::Val{:y})::Vector{T} where {T <: SSDFloat}
    return geom_round.(T[-hp.rOuter, -sqrt(3)/2*hp.rOuter, -hp.rInner, -sqrt(3)/2*hp.rInner, 0., sqrt(3)/2*hp.rInner, hp.rInner, sqrt(3)/2*hp.rOuter, hp.rOuter] .+ hp.translate.y)
end

#and a sample function to paint the primitive on the grid (necessary if the object is small)
function sample(hp::HexagonalPrism{T}, stepsize::Vector{T})  where {T <: SSDFloat}
    samples = CartesianPoint{T}[]
    for x in -hp.rOuter : stepsize[1] : hp.rOuter
        for y in -hp.rOuter : stepsize[2]: hp.rOuter
            if CartesianPoint{T}(x, y, 0) in hp
                for z in -hp.h/2 : stepsize[3] : hp.h/2
                    push!(samples, CartesianPoint{T}(x, y, z))
                end
            end
        end
    end
    ismissing(hp.translate) ? nothing : samples = map(x -> x + hp.translate, samples)
    return samples
end
