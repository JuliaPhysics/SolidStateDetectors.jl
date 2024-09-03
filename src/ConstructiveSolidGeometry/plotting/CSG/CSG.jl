@inline primitives(vp::AbstractVolumePrimitive) = (vp,)
@inline function primitives(csg::AbstractConstructiveGeometry)
    (
        primitives(csg.a)...,
        primitives(csg.b)...
    )
end

function _get_csg_plot_objects(csg::AbstractConstructiveGeometry{T}, primitive, n_samples, CSG_scale, slice_plane, st) where {T}
    if st == :samplesurface
        CSG_scale = ismissing(CSG_scale) ? get_scale(csg) : CSG_scale
        spacing = T(CSG_scale/n_samples)
        filter(p -> in(p,csg), sample(primitive, spacing))
    elseif st == :slice
        axis, slice  = slice_plane
        CSG_scale = ismissing(CSG_scale) ? get_scale(csg) : CSG_scale
        pointtype, axisrot = CartesianPoint, axis
        spacing = T(CSG_scale/n_samples)
        if axis == :φ
            pointtype, axisrot, slice, spacing = CylindricalPoint, :y, 0, spacing*T(0.5)
        end
        samples = filter(pt -> abs(getproperty(pt, axisrot) - T(to_internal_units(slice))) < spacing/2, filter(p -> in(p,csg), sample(primitive, spacing)))
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
            _get_csg_plot_objects(csg::AbstractConstructiveGeometry{T}, primitive, n_samples, CSG_scale, slice_plane, st)
        else
            internal_length_unit*u, internal_length_unit*v
        end
    else
        primitive
    end
end

@recipe function f(csg::AbstractConstructiveGeometry{T}; n_samples = 40, CSG_scale = missing, slice_plane = (:y, 0)) where {T}
    seriestype --> :csg
    st = plotattributes[:seriestype]
    unitformat --> :slash
    axis, slice  = slice_plane
    pointtype = axis == :φ ? CylindricalPoint : CartesianPoint
    if st == :slice && axis == :φ
        csg = transform(csg, (rotation = T.(RotZ(-to_internal_units(slice))), translation = CartesianVector{T}()))
    end
    ps = primitives(csg)
    
    @series begin
        label --> "CSG"
        if st == :slice
            proj = filter(x -> x != axis, fieldnames(pointtype))
            xunit --> internal_length_unit
            yunit --> internal_length_unit
            xguide --> string(proj[1])
            yguide --> string(proj[2])
        else
            if !isClosedPrimitive(ps[1])
                fillcolor := :white
                fillalpha --> 0.2
                linewidth := 1
            end
        end
        _get_csg_plot_objects(csg, ps[1], n_samples, CSG_scale, slice_plane, st)
    end
    for i in 2:length(ps)
        @series begin
            label := ""
            if st == :slice
                proj = filter(x -> x != axis, fieldnames(pointtype))
                xunit --> internal_length_unit
                yunit --> internal_length_unit
                xguide --> string(proj[1])
                yguide --> string(proj[2])
            else
                if !isClosedPrimitive(ps[i])
                    fillcolor := :white
                    fillalpha --> 0.2
                    linewidth := 1
                end
            end
            _get_csg_plot_objects(csg, ps[i], n_samples, CSG_scale, slice_plane, st)
        end
    end
end