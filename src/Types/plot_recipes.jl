@recipe function f( pts::PointTypes;
                    r=missing,
                    φ=missing,
                    z=missing )
    T = eltype(pts.r)
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

    maximum_voltage::T = maximum(pts.pointtypes)
    minimum_voltage::T = minimum(pts.pointtypes)

    length_r::T = pts.r[end] - pts.r[1]
    length_φ::T = pts.φ[end] - pts.φ[1]
    length_z::T = pts.z[end] - pts.z[1]

    clims --> (minimum_voltage, maximum_voltage)

    if cross_section == :φ
        label --> ""
        xlabel --> L"$r$ / m"
        ylabel --> L"$z$ / m"

        aspect_ratio --> 1
        plt_ratio = length_r / length_z
        size --> (1200 * plt_ratio, 1000)

        iφ = searchsortedfirst(pts.φ, value )
        if iφ > length(pts.φ) iφ = length(pts.φ) end

        title --> L"\varphi" * " = $(round(pts.φ[iφ], sigdigits=3)) rad ($(round(rad2deg(pts.φ[iφ]), sigdigits=3))" * L"\degree" * ")"

        pts.r, pts.z, pts.pointtypes[:,iφ,:]'
    elseif cross_section == :z
        proj --> :polar
        aspect_ratio --> 1
        ylims --> (0., pts.r[end]*1.05)

        iz = searchsortedfirst(pts.z, value )
        if iz > length(pts.z) iz = length(pts.z) end

        size --> (900, 900)

        title --> L"z" * " = $(round(1e3 * pts.z[iz], sigdigits=3)) mm"

        if length(pts.φ) < 2
            error(" Grid has only 1 tick in φ. Cannot be displayed in φ-r plane.\n Extend it via `SolidStateDetectors.extent_2D_grid_to_3D!(g::CylindricalGrid, n::Int)` to n ticks in φ.")
        end          

        pts.φ, pts.r, pts.pointtypes[:, :, iz]
    elseif cross_section == :r
        label --> ""
        xlabel --> L"\varphi ⋅ r" * " / mm"
        ylabel --> L"z" * " / mm"

        aspect_ratio --> 1

        ir = searchsortedfirst(pts.r, value )
        if ir > length(pts.r) ir = length(pts.r) end

        title --> L"r" * " = $(round(1e3 * pts.r[ir], sigdigits=3)) mm"

        plt_ratio = length_φ * pts.r[ir] / length_z
        size --> (1200 * plt_ratio, 1000)

        pts.φ * pts.r[ir], pts.z, pts.pointtypes[ir, :, :]'
    end
end
