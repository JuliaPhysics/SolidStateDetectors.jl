function sample(p::AbstractVolumePrimitive{T}, spacing::T)::Vector{CartesianPoint{T}} where {T}
    fs = surfaces(p)
    vs = Vector{CartesianPoint{T}}()
    for s in fs
        append!(vs, sample(s, spacing))
    end
    vs
end

function get_crosssection_and_value(x::Union{Real, LengthQuantity, Missing}, y::Union{Real, LengthQuantity, Missing}, z::Union{Real, LengthQuantity, Missing}, φ::Union{Real, AngleQuantity, Missing})::Tuple{Symbol, Real}
    crosssections = Dict(:x => x, :y => y, :z => z, :φ => φ)
    idx = findfirst(cs -> !ismissing(cs), crosssections)
    if isnothing(idx)
        :y, 0
    else
        return idx, to_internal_units(crosssections[idx])
    end
end

function is_coplanar(s::AbstractSurfacePrimitive{T}, plane::Tuple{Symbol,T}) where {T}
    axis, slice = plane
    crosssections = Dict(   :x => CartesianVector{T}(1,0,0), 
                            :y => CartesianVector{T}(0,1,0), 
                            :z => CartesianVector{T}(0,0,1),
                            :φ => CartesianVector{T}(cos(π/2+slice),sin(π/2+slice),0)
                   )
    pointtype = axis == :φ ? CylindricalPoint : CartesianPoint
    point = pointtype(vertices(s,1)[1])
    abs(abs(dot(normal(s), crosssections[axis])) - 1) < csg_default_tol(T) && abs(getproperty(point, axis) - slice) < csg_default_tol(T)
end

@recipe function f(p::AbstractPrimitive{T}, axis::Symbol, slice::T, st::Symbol; n_samples = 40, CSG_scale = missing, projection = :none, full_det = false, linewidth = :auto) where {T}
    if st == :samplesurface
        CSG_scale = ismissing(CSG_scale) ? extremum(p) : CSG_scale
        spacing = T(CSG_scale/n_samples)
        sample(p, spacing)
    elseif st == :slice
        isgr = occursin("GRBackend", string(typeof(plotattributes[:plot_object].backend)))
        pointtype, axisrot = CartesianPoint, axis
        CSG_scale = ismissing(CSG_scale) ? extremum(p) : CSG_scale
        spacing = T(CSG_scale/n_samples)
        if axis == :φ
            p = rotate(p, T.(RotZ(-to_internal_units(slice))))
            pointtype, axisrot, slice, spacing = CylindricalPoint, :y, 0, spacing*T(0.5)
        end
        samples = filter(pt -> abs(getproperty(pt, axisrot) - T(to_internal_units(slice))) < spacing/2, sample(p, spacing))
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
        [surfaces(p)...]
    end
end

@recipe function f(p::AbstractVolumePrimitive{T}; x = missing, y = missing, z = missing, φ = missing) where {T}
    linecolor --> :black
    seriestype --> :csg
    st = plotattributes[:seriestype]
    unitformat --> :slash
    axis, slice  = get_crosssection_and_value(x,y,z,φ)
    if st == :slice
        surfs = filter(s -> typeof(s) <: AbstractPlanarSurfacePrimitive{T}, surfaces(p))
        if any(s -> is_coplanar(s, (axis, T(slice))), surfs)
            n_samples --> 100
        else
            n_samples --> 200
        end
    else
        n_samples --> 40
    end
    @series begin
        label --> "$(nameof(typeof(p)))"
        p, axis, T(slice), st
    end
end