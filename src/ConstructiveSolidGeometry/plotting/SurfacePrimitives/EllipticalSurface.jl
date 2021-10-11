function mesh(es::EllipticalSurface{T}; n = 30)::Mesh{T} where {T}

    rMin::T, rMax::T = get_r_limits(es)
    φMin::T, φMax::T = get_φ_limits(es)
    f = (φMax - φMin)/(2π)
    n = Int(ceil(n*f))
    φ = range(φMin, φMax, length = n+1)
    r = [rMin, rMax]
    if rMin == 0
        x = append!([0.0],rMax*cos.(φ))
        y = append!([0.0],rMax*sin.(φ))
        z = zeros(n+2)
        connections = [[1,i,i+1] for i in 2:n+1]
    else
        x = append!(rMin*cos.(φ), rMax*cos.(φ))
        y = append!(rMin*sin.(φ), rMax*sin.(φ))
        z = zeros(2(n+1))
        connections = [[i,i+1,i+n+2,i+n+1] for i in 1:n]
    end

    es.rotation*Mesh{T}(x,y,z,connections) + es.origin
end

@recipe function f(es::EllipticalSurface; n = 40)
    seriestype --> :mesh3d
    if haskey(plotattributes, :seriestype) && plotattributes[:seriestype] == :mesh3d
        @series begin
            label --> "Elliptical Surface"
            mesh(es, n = n)
        end
    else
        ls = lines(es)
        linecolor --> :black
        @series begin
            label --> "Elliptical Surface"
            n := n
            ls[1]
        end
        if length(ls) > 1 
            for i in 2:length(ls)
                @series begin 
                    label := nothing
                    n := n
                    ls[i]
                end
            end
        end
    end
    if haskey(plotattributes, :show_normal) && plotattributes[:show_normal]
        @series begin
            label := nothing
            seriestype := :vector
            _plt_get_start_point_for_normal(es), normalize(Plane(es).normal) * sqrt(_plt_area(es)) / 20 
        end
    end
end

_plt_area(es::CircularArea{T}) where {T} = π*es.r^2
_plt_area(es::PartialCircularArea{T}) where {T} = π*es.r^2 / ((es.φ[2]-es.φ[1])/2π)
_plt_area(es::Annulus{T}) where {T} = π*(es.r[2]^2 - es.r[1]^2)
_plt_area(es::PartialAnnulus{T}) where {T} = π*(es.r[2]^2 - es.r[1]^2) / ((es.φ[2]-es.φ[1])/2π)

_plt_get_start_point_for_normal(es::CircularArea{T}) where {T} = es.origin
function _plt_get_start_point_for_normal(es::PartialCircularArea{T}) where {T}
    cyl = CylindricalPoint{T}(es.r/2, (es.φ[2]+es.φ[1])/2, zero(T))
    _transform_into_global_coordinate_system(CartesianPoint(cyl), es)
end
function _plt_get_start_point_for_normal(es::Annulus{T}) where {T}
    cyl = CylindricalPoint{T}((es.r[2]+es.r[1])/2, zero(T), zero(T))
    _transform_into_global_coordinate_system(CartesianPoint(cyl), es)
end
function _plt_get_start_point_for_normal(es::PartialAnnulus{T}) where {T}
    cyl = CylindricalPoint{T}((es.r[2]+es.r[1])/2, (es.φ[2]+es.φ[1])/2, zero(T))
    _transform_into_global_coordinate_system(CartesianPoint(cyl), es)
end
