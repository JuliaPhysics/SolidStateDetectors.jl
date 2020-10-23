
function LineSegments(t::Tube{T})::Vector{AbstractLine{T,3,:cartesian}} where {T <: SSDFloat}
    ls = AbstractLine{T, 3, :cartesian}[]
    translate::CartesianVector{T} = ismissing(t.translate) ? CartesianVector{T}([0, 0, 0]) : t.translate
    for r in (t.r_interval.left == 0 ? [t.r_interval.right] : [t.r_interval.left, t.r_interval.right])
        for z in [t.z_interval.left, t.z_interval.right]
            push!(ls, PartialCircle(r, t.φ_interval.left, t.φ_interval.right, translate + CartesianVector{T}([0, 0, z])))
        end
    end
    for r in [t.r_interval.left, t.r_interval.right]
        if r != 0
            for φ in ((t.φ_interval.right - t.φ_interval.left ≈ 2π) ? [t.φ_interval.left] : [t.φ_interval.left, t.φ_interval.right])
                push!(ls, LineSegment(
                    CartesianPoint{T}(r * cos(φ), r * sin(φ), t.z_interval.left) + translate,
                    CartesianPoint{T}(r * cos(φ), r * sin(φ), t.z_interval.right) + translate))
            end
        end
    end
    for φ in ((t.φ_interval.right - t.φ_interval.left ≈ 2π) ? [] : [t.φ_interval.left, t.φ_interval.right])
        for z in [t.z_interval.left, t.z_interval.right]
            push!(ls, LineSegment(
                CartesianPoint{T}(t.r_interval.left * cos(φ), t.r_interval.left * sin(φ), z) + translate,
                CartesianPoint{T}(t.r_interval.right * cos(φ), t.r_interval.right * sin(φ), z) + translate))
        end
    end
    return ls
end

@recipe function f(t::Tube{T}; n = 30, seriescolor = :green) where {T}
    linewidth --> 2
    n --> n
    @series begin
        seriescolor --> seriescolor
        label --> "Tube"
        []
    end
    label := ""
    seriescolor := seriescolor
    LineSegments(t)
end


function LineSegments(c::Cone{T})::Vector{AbstractLine{T, 3, :cartesian}} where {T <: SSDFloat}
    ls = AbstractLine{T, 3, :cartesian}[]
    translate::CartesianVector{T} = ismissing(c.translate) ? CartesianVector{T}([0, 0, 0]) : c.translate
    for r in [c.rStart1, c.rStop1]
        push!(ls, PartialCircle(r, c.φStart, c.φStop, translate + CartesianVector{T}([0, 0, c.zStart]))) 
    end
    for r in [c.rStart2, c.rStop2]
        push!(ls, PartialCircle(r, c.φStart, c.φStop, translate + CartesianVector{T}([0, 0, c.zStop])))

    end
    for φ in ((c.φStop - c.φStart ≈ 2π) ? [c.φStart] : [c.φStart, c.φStop])
        push!(ls, LineSegment(
            CartesianPoint{T}(c.rStart1 * sin(φ), c.rStart1 * cos(φ), c.zStart) + translate,
            CartesianPoint{T}(c.rStart2 * sin(φ), c.rStart2 * cos(φ), c.zStop) + translate))
        push!(ls, LineSegment(
            CartesianPoint{T}(c.rStop1 * sin(φ), c.rStop1 * cos(φ), c.zStart) + translate,
            CartesianPoint{T}(c.rStop2 * sin(φ), c.rStop2 * cos(φ), c.zStop) + translate))
    end
    for φ in ((c.φStop - c.φStart ≈ 2π) ? [c.φStart] : [c.φStart, c.φStop])
        push!(ls, LineSegment(
            CartesianPoint{T}(c.rStart1 * sin(φ),  c.rStart1 * cos(φ),  c.zStart) + translate,
            CartesianPoint{T}(c.rStop1 * sin(φ),   c.rStop1 * cos(φ), c.zStart) + translate))
        push!(ls, LineSegment(
            CartesianPoint{T}(c.rStart2 * sin(φ),  c.rStart2 * cos(φ),  c.zStop) + translate,
            CartesianPoint{T}(c.rStop2 * sin(φ),   c.rStop2 * cos(φ), c.zStop) + translate))
    end
    return ls
end

@recipe function f(c::Cone{T}; n = 30, seriescolor = :orange) where {T}
    linewidth --> 2
    n --> n
    @series begin
        seriescolor --> seriescolor
        label --> "Cone"
        []
    end
    seriescolor := seriescolor
    label := ""
    LineSegments(c)
end