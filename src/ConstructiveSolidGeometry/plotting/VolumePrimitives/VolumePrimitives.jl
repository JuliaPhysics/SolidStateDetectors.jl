function sample(p::AbstractVolumePrimitive{T}, spacing::T)::Vector{CartesianPoint{T}} where {T}
    fs = surfaces(p)
    vs = Vector{CartesianPoint{T}}()
    for s in fs
        append!(vs, sample(s, spacing))
    end
    vs
end

function _get_volume_primitive_plot_objects(p::AbstractVolumePrimitive{T}, n_samples, slice_plane, st) where {T}
    if st == :samplesurface
        spacing = T(extremum(p)/n_samples)
        sample(p, spacing)
    elseif st == :slice
        axis, slice  = slice_plane
        pointtype, axisrot = CartesianPoint, axis
        spacing = T(extremum(p)/n_samples)
        if axis == :φ
            p = rotate(p, T.(RotZ(-to_internal_units(slice))))
            pointtype, axisrot, slice, spacing = CylindricalPoint, :y, 0, spacing*T(0.5)
        end
        samples = filter(pt -> abs(getproperty(pt, axisrot) - T(to_internal_units(slice))) < spacing/2, sample(p, spacing))
        proj = fieldnames(CartesianPoint)[findall(x -> x != axis, fieldnames(pointtype))]
        u = getproperty.(samples, proj[1])
        v = getproperty.(samples, proj[2])
        if axis == :φ
            idx = findall(0 .≤ u)
            u = u[idx]
            v = v[idx]
        end
        if length(u) > 10000
            n_samples = Int(floor(n_samples*0.5))
            _get_volume_primitive_plot_objects(p, n_samples, slice_plane, st)
        else
            internal_length_unit*u, internal_length_unit*v
        end
    else
        [surfaces(p)...]
    end
end

@recipe function f(p::AbstractVolumePrimitive{T}; n_samples = 40, slice_plane = (:y, 0)) where {T}
    linecolor --> :black
    seriestype --> :csg
    st = plotattributes[:seriestype]
    axis, slice  = slice_plane
    pointtype = axis == :φ ? CylindricalPoint : CartesianPoint
    @series begin
        label --> "$(nameof(typeof(p)))"
        if st == :slice
            proj = filter(x -> x != axis, fieldnames(pointtype))
            xguide --> string(proj[1])
            yguide --> string(proj[2])
            xunit --> internal_length_unit
            yunit --> internal_length_unit
            unitformat --> :slash
        end
        _get_volume_primitive_plot_objects(p, n_samples, slice_plane, st)
    end
end