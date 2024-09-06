@inline primitives(vp::AbstractVolumePrimitive) = (vp,)
@inline function primitives(csg::AbstractConstructiveGeometry)
    (
        primitives(csg.a)...,
        primitives(csg.b)...
    )
end

@recipe function f(csg::AbstractConstructiveGeometry{T}, primitive::AbstractPrimitive{T}, axis::Symbol, slice::T, st::Symbol; n_samples = 40, CSG_scale = missing, projection = :none, full_det = false, linewidth = :auto) where {T}
    if st == :samplesurface
        CSG_scale = ismissing(CSG_scale) ? get_scale(csg) : CSG_scale
        spacing = T(CSG_scale/n_samples)
        filter(p -> in(p,csg), sample(primitive, spacing))
    elseif st == :slice
        isgr = occursin("GRBackend", string(typeof(plotattributes[:plot_object].backend)))
        CSG_scale = ismissing(CSG_scale) ? get_scale(csg) : CSG_scale
        pointtype, axisrot = CartesianPoint, axis
        spacing = T(CSG_scale/n_samples)
        if axis == :φ
            pointtype, axisrot, slice, spacing = CylindricalPoint, :y, 0, spacing*T(0.5)
        end
        samples = filter(pt -> in(pt, csg) && abs(getproperty(pt, axisrot) - slice) < spacing/2, sample(primitive, spacing))
        proj = fieldnames(CartesianPoint)[findall(x -> x != axis, fieldnames(pointtype))]
        u = getproperty.(samples, proj[1])
        v = getproperty.(samples, proj[2])
        if axis == :φ && !full_det
            idx = findall(0 .≤ u)
            u = u[idx]
            v = v[idx]
        end
        proj = filter(x -> x != axis, fieldnames(pointtype))
        yunit --> internal_length_unit
        xguide --> string(proj[1])
        yguide --> string(proj[2])
        if linewidth != :auto
            markersize := linewidth
        end
        if projection == :polar && axis == :z
            xunit --> internal_angle_unit
            uconvert.(internal_angle_unit, atan.(v,u)), internal_length_unit*sqrt.(u.^2 + v.^2)
        else
            if isgr aspect_ratio --> 1.0 end
            xunit --> internal_length_unit
            internal_length_unit*u, internal_length_unit*v
        end
    else
        if !isClosedPrimitive(primitive)
            fillcolor := :white
            fillalpha --> 0.2
            linewidth := 1
        end
        primitive
    end
end

@recipe function f(csg::AbstractConstructiveGeometry{T}; x = missing, y = missing, z = missing, φ = missing) where {T}
    seriestype --> :csg
    st = plotattributes[:seriestype]
    unitformat --> :slash
    axis, slice  = get_crosssection_and_value(x,y,z,φ)
    if st == :slice
        surfs = filter(s -> typeof(s) <: AbstractPlanarSurfacePrimitive{T}, surfaces(csg))
        if any(s -> is_coplanar(s, (axis, T(slice))), surfs)
            n_samples --> 100
        else
            n_samples --> 200
        end
        if axis == :φ
            csg = transform(csg, (rotation = T.(RotZ(-to_internal_units(slice))), translation = CartesianVector{T}()))
        end
    else
        n_samples --> 40
    end
    ps = primitives(csg)
    
    @series begin
        label --> "CSG"
        csg, ps[1], axis, T(slice), st
    end
    for i in 2:length(ps)
        @series begin
            label := ""
            csg, ps[i], axis, T(slice), st
        end
    end
end