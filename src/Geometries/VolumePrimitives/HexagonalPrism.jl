struct HexagonalPrism{T} <: AbstractVolumePrimitive{T, 3} ## Only upright hexagons at the moment
    org::CartesianPoint{T} #origin at middle of hexagonal prism
    a::T #side length/ maximal radius
    h::T #total height
end

# is a point inside the hexagonal prism?
function in(pt::CartesianPoint{T}, hp::HexagonalPrism{T})::Bool where {T}
    shift::CartesianPoint{T} = pt - hp.org #shift pt to prism's frame
    r = cos(deg2rad(30))*hp.a #minimal radius
    #convert to cylindrical to break into 6 segements of a hexagon
    shift2 = CylindricalPoint(shift)
    θ = rad2deg(shift2.φ)
    ρ = shift2.r
    return  abs(shift.z) <= hp.h/2 &&
            ρ <= hp.a && #outside maximal radius?
            if ρ > r     #if outside minimal radius we need to carefully check each segment
                if (θ > 60 && θ < 120)
                    shift.y < r
                elseif (θ > 240 && θ < 300)
                    shift.y > -r
                #for the rest we rotate so we only need to compare y
                elseif (θ > 0 && θ < 60)
                    y_rot = shift.x*sin(deg2rad(60)) + shift.y*cos(deg2rad(60))
                    y_rot < r
                elseif (θ > 120 && θ < 180)
                    y_rot = shift.x*sin(deg2rad(-60)) + shift.y*cos(deg2rad(-60))
                    y_rot < r
                elseif (θ > 180 && θ < 240)
                    y_rot = shift.x*sin(deg2rad(60)) + shift.y*cos(deg2rad(60))
                    y_rot > -r
                elseif (θ > 300 && θ < 360)
                    y_rot = shift.x*sin(deg2rad(-60)) + shift.y*cos(deg2rad(-60))
                    y_rot > -r
                else
                    true
                end
            else
                true
            end
    end

function in(pt::CylindricalPoint{T}, hp::HexagonalPrism{T})::Bool where {T}
    return in(CartesianPoint(pt), hp)
end


# You also have to implement the function to obtain the primitive from a config file (so an dic)
# You also should provide a example config file containing this new primitive
function HexagonalPrism{T}(dict::Union{Dict{Any, Any}, Dict{String, Any}}, inputunit_dict::Dict{String,Unitful.Units})::HexagonalPrism{T} where {T <: SSDFloat}
    # ... parse values from dict to NewPrimitive{T}(...)
    if haskey(dict, "translate")
        org::CartesianPoint{T} = CartesianPoint{T}(
            haskey(dict["translate"],"x") ? geom_round(ustrip(uconvert(u"m", T(dict["translate"]["x"]) * inputunit_dict["length"] ))) : 0,
            haskey(dict["translate"],"y") ? geom_round(ustrip(uconvert(u"m", T(dict["translate"]["y"]) * inputunit_dict["length"] ))) : 0,
            haskey(dict["translate"],"z") ? geom_round(ustrip(uconvert(u"m", T(dict["translate"]["z"]) * inputunit_dict["length"] ))) : 0)
    else
        org = CartesianPoint{T}(0,0,0)
    end
    a::T = ustrip(uconvert(u"m", dict["a"] * inputunit_dict["length"]))
    h::T = ustrip(uconvert(u"m", dict["h"] * inputunit_dict["length"]))
    #φ::T = ustrip(uconvert(u"m", dict["φ"] * inputunit_dict["length"]))
    return HexagonalPrism{T}(org, a, h)
end
function Geometry(T::DataType, t::Val{:HexagonalPrism}, dict::Dict{Any, Any}, inputunit_dict::Dict{String,Unitful.Units})
    return HexagonalPrism{T}(dict, inputunit_dict)
end

# add a (+) method to shift the primitive
function (+)(hp::HexagonalPrism{T}, translate::Union{CartesianVector{T}, Missing})::HexagonalPrism{T} where {T <: SSDFloat}
    return ismissing(translate) ? hp : HexagonalPrism{T}( hp.org + translate, hp.a, hp.h)
end

# Also a plot recipe for this new primitive should be provided:
@recipe function f(hp::HexagonalPrism{T}) where {T <: SSDFloat}
    label --> "HexagonalPrism"
    @series begin
        pts_top = []
        pts_bottom = []

        #find all vertices, this loop has been tested and works
        for φ in [0,deg2rad(60), deg2rad(120), deg2rad(180), deg2rad(240), deg2rad(300)]
            pts_top = hp.org + CartesianPoint(hp.a * cos(φ), hp.a * sin(φ), hp.h/2)
            push!(pts_top, pt_top)
            pts_bottom = hp.org + CartesianPoint(hp.a * cos(φ), hp.a * sin(φ), -hp.h/2)
            push!(pts_bottom, pt_bottom)
        end

        #create Linesegments connecting the vertices
        lines = LineSegment{T, 3, :cartesian}[]
        for i in 1:length(pts_top)-1
            push!(lines, LineSegment(pts_top[i+1], pts_top[i])) #top hexagon
            push!(lines, LineSegment(pts_bottom[i+1], pts_bottom[i])) #bottom hexagon
            push!(lines, LineSegment(pts_bottom[i], pts_top[i])) #lines connecting hexagons (missing one line)
        lines
    end
end

# For proper grid creation we also need the function get_important_points:
function get_important_points(hp::HexagonalPrism{T}, ::Val{:r})::Vector{T} where {T <: SSDFloat}
    r = hp.a > (hp.h/2) ? hp.a : hp.h/2
    v::Vector{T} = geom_round.(T[ hp.org.x - r, hp.org.x, hp.org.x + r,
                                  hp.org.y - r, hp.org.y, hp.org.y + r ])
    return findall(r -> r >= 0 , v)
end
function get_important_points(hp::HexagonalPrism{T}, ::Val{:φ})::Vector{T} where {T <: SSDFloat}
    return T[]
end
function get_important_points(hp::HexagonalPrism{T}, ::Val{:z})::Vector{T} where {T <: SSDFloat}
    return geom_round.(T[hp.org.z+hp.h/2, hp.org.z-hp.h/2])
end
function get_important_points(hp::HexagonalPrism{T}, ::Val{:x})::Vector{T} where {T <: SSDFloat}
    return geom_round.(T[hp.org.x-hp.a,hp.org.x+hp.a])
end
function get_important_points(hp::HexagonalPrism{T}, ::Val{:y})::Vector{T} where {T <: SSDFloat}
    return geom_round.(T[hp.org.y-hp.a, hp.org.y+hp.a])
end

and a sample function to paint the primitive on the grid (necessary if the object is small)
function sample(hp::HexagonalPrism{T}, stepsize::Vector{T})  where {T <: SSDFloat}
    samples::Vector{CartesianPoint{T}} = CartesianPoint{T}[]
    φarr::Vector{T} = geom_round.(T[0,deg2rad(60), deg2rad(120), deg2rad(180), deg2rad(240), deg2rad(300)])
    r = cos(deg2rad(30))*hp.a
    zarr::Vector{T} = get_important_points(hp, Val{:z}())
    for φ in φarr
        for z in zarr[1] : zarr[2]
            pt1::CylindricalPoint{T} = CylindricalPoint{T}(cos(φ)*hp.a+hp.org.x, sin(φ)*hp.a+hp.org.y, z)
            pt2::CylindricalPoint{T} = CylindricalPoint{T}(cos(φ+deg2rad(30))*r+hp.org.x, sin(φ+deg2rad(30))*r+hp.org.y, z)
            push!(samples, pt1)
            push!(samples, pt2)
        end
    end
    return samples
end
