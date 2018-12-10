@recipe function f( g::CylindricalGrid;
                    r=missing,
                    φ=missing,
                    z=missing,
                    contours_equal_potential=false )
    T = eltype(g.r)
    seriescolor --> :viridis
    st --> :heatmap

    cross_section::Symbol, value::T = if ismissing(φ) && ismissing(r) && ismissing(z)
        :φ, 0
    elseif !ismissing(φ) && ismissing(r) && ismissing(z)
        :φ, φ
    elseif ismissing(φ) && !ismissing(r) && ismissing(z)
        :r, r
    elseif ismissing(φ) && ismissing(r) && !ismissing(z)
        :z, z
    else
        error(ArgumentError, ": Only one of the keywords `r, φ, z` is allowed.")
    end

    maximum_voltage::T = maximum(g.potential)
    minimum_voltage::T = minimum(g.potential)

    length_r::T = g.r[end] - g.r[1]
    length_φ::T = g.φ[end] - g.φ[1]
    length_z::T = g.z[end] - g.z[1]

    clims --> (minimum_voltage, maximum_voltage)

    if cross_section == :φ
        label --> ""
        xlabel --> L"$r$ / m"
        ylabel --> L"$z$ / m"

        aspect_ratio --> 1
        plt_ratio = length_r / length_z
        size --> (700 * plt_ratio + 100, 700)

        iφ = searchsortedfirst(g.φ, value )
        if iφ > length(g.φ) iφ = length(g.φ) end
        title --> L"\varphi" * " = $(round(g.φ[iφ], sigdigits=3)) rad ($(round(rad2deg(g.φ[iφ]), sigdigits=3))" * L"\degree" * ")"

        @series begin
            g.r, g.z, g.potential[:,iφ,:]'
        end
        if contours_equal_potential
            @series begin
                seriescolor := :thermal
                st := :contours
                g.r, g.z, g.potential[:,iφ,:]'
            end
        end
    elseif cross_section == :z
        proj --> :polar
        aspect_ratio --> 1
        ylims --> (0., g.r[end]*1.05)

        iz = searchsortedfirst(g.z, value )
        if iz > length(g.z) iz = length(g.z) end

        title --> L"z" * " = $(round(1e3 * g.z[iz], sigdigits=3)) mm"
        size --> (800, 800)

        if length(g.φ) < 2
            error(" Grid has only 1 tick in φ. Cannot be displayed in φ-r plane.\n Extend it via `SolidStateDetectors.extent_2D_grid_to_3D!(g::CylindricalGrid, n::Int)` to n ticks in φ.")
        end        
        @series begin
            g.φ, g.r, g.potential[:, :, iz]
        end
        if contours_equal_potential
            @series begin
                seriescolor := :thermal
                st := :contours
                vcat(g.φ, T[2π]), g.r, cat(g.potential[:, :, iz], g.potential[:, 1, iz], dims=2)
            end
        end
    elseif cross_section == :r
        label --> ""
        xlabel --> L"\varphi ⋅ r" * " / mm"
        ylabel --> L"z" * " / mm"

        aspect_ratio --> 1

        ir = searchsortedfirst(g.r, value )
        if ir > length(g.r) ir = length(g.r) end

        title --> L"r" * " = $(round(1e3 * g.r[ir], sigdigits=3)) mm"

        plt_ratio = length_φ * g.r[ir] / length_z
        if plt_ratio == 0 plt_ratio = 1 end
        size --> (600 * plt_ratio + 100, 600)
        @series begin
            g.φ * g.r[ir], g.z, g.potential[ir, :, :]'
        end
        if contours_equal_potential
            @series begin
                levels --> 20
                seriescolor := :thermal
                st := :contours
                g.φ * g.r[ir], g.z, g.potential[ir, :, :]'
            end
        end
    end
end

@userplot GridDensityDisplay
@recipe function f(gdd::GridDensityDisplay)
    if !(   typeof(gdd.args[1]) == CylindricalGrid{Float32} || typeof(gdd.args[1]) == CylindricalGrid{Float64}
            || typeof(gdd.args[1]) == CylindricalGrid{Float16} )
        error("First input argument should be of type `SolidStateDetectors.DetectorGrid.CylindricalGrid`.")
    end
    cg = gdd.args[1]

    size --> (1200, 400)
    layout := (1, 3)
    legend --> false

    @series begin
        subplot := 1
        xlabel := "r / m"
        seriestype --> :stephist
        title --> "$(length(cg.r)) points in r"
        cg.r
    end

    @series begin
        subplot := 2
        xlabel := "φ / rad"
        seriestype --> :stephist
        title --> "$(length(cg.φ)) points in φ"
        cg.φ
    end

    @series begin
        subplot := 3
        xlabel := "z / m"
        seriestype --> :stephist
        title --> "$(length(cg.z)) points in z"
        cg.z
    end
end
