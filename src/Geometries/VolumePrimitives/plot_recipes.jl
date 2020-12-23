
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
            CartesianPoint{T}(c.rStart1 * cos(φ), c.rStart1 * sin(φ), c.zStart) + translate,
            CartesianPoint{T}(c.rStart2 * cos(φ), c.rStart2 * sin(φ), c.zStop) + translate))
        push!(ls, LineSegment(
            CartesianPoint{T}(c.rStop1 * cos(φ), c.rStop1 * sin(φ), c.zStart) + translate,
            CartesianPoint{T}(c.rStop2 * cos(φ), c.rStop2 * sin(φ), c.zStop) + translate))
    end
    for φ in ((c.φStop - c.φStart ≈ 2π) ? [c.φStart] : [c.φStart, c.φStop])
        push!(ls, LineSegment(
            CartesianPoint{T}(c.rStart1 * cos(φ),  c.rStart1 * sin(φ),  c.zStart) + translate,
            CartesianPoint{T}(c.rStop1 * cos(φ),   c.rStop1 * sin(φ), c.zStart) + translate))
        push!(ls, LineSegment(
            CartesianPoint{T}(c.rStart2 * cos(φ),  c.rStart2 * sin(φ),  c.zStop) + translate,
            CartesianPoint{T}(c.rStop2 * cos(φ),   c.rStop2 * sin(φ), c.zStop) + translate))
    end
    return ls
end



function LineSegments(t::Torus{T})::Vector{AbstractLine{T,3,:cartesian}} where {T <: SSDFloat}
    ls = AbstractLine{T, 3, :cartesian}[]
    translate::CartesianVector{T} = ismissing(t.translate) ? CartesianVector{T}([0, 0, 0]) : t.translate
    for r_tube in (t.r_tube_interval.left == 0 ? [t.r_tube_interval.right] : [t.r_tube_interval.left, t.r_tube_interval.right])
        θ_circles = (t.θ_interval.right - t.θ_interval.left) ≈ 2π ? [t.θ_interval.left] : [t.θ_interval.left, t.θ_interval.right]
        π in OpenInterval(t.θ_interval) ? push!(θ_circles, π) : nothing
        for θ in θ_circles
            r = r_tube*cos(θ) + t.r_torus
            z = r_tube*sin(θ)
            push!(ls, PartialCircle(r, t.φ_interval.left, t.φ_interval.right, translate + CartesianVector{T}([0, 0, z])))
        end
    end
    for φ in (t.φ_interval.right - t.φ_interval.left ≈ 2π ? [t.φ_interval.left] : [t.φ_interval.left, t.φ_interval.right])
        for r_tube in (t.r_tube_interval.left == 0 ? [t.r_tube_interval.right] : [t.r_tube_interval.left, t.r_tube_interval.right])
            r = r_tube
            push!(ls, PartialCircle(r, t.θ_interval.left, t.θ_interval.right, translate + CartesianVector{T}([t.r_torus*cos(φ), t.r_torus*sin(φ), 0]), RotZ{T}(φ)*RotX{T}(π/2)))
        end
        for θ in (t.θ_interval.right - t.θ_interval.left ≈ 2π ? [] : [t.θ_interval.left, t.θ_interval.right])
            push!(ls, LineSegment(
                RotZ{T}(φ)*RotY{T}(-θ)*CartesianPoint{T}(t.r_tube_interval.left, 0, 0) + CartesianVector{T}([t.r_torus*cos(φ), t.r_torus*sin(φ), 0]) + translate,
                RotZ{T}(φ)*RotY{T}(-θ)*CartesianPoint{T}(t.r_tube_interval.right, 0, 0) + CartesianVector{T}([t.r_torus*cos(φ), t.r_torus*sin(φ), 0]) + translate))
        end
    end
    return ls
end

