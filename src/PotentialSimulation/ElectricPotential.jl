struct ElectricPotential{T, N, S} <: AbstractArray{T, N}
    data::Array{T, N}
    grid::Grid{T, N, S}
end

@inline size(ep::ElectricPotential{T, N, S}) where {T, N, S} = size(ep.data)
@inline length(ep::ElectricPotential{T, N, S}) where {T, N, S} = length(ep.data)
@inline getindex(ep::ElectricPotential{T, N, S}, I::Vararg{Int, N}) where {T, N, S} = getindex(ep.data, I...)
@inline getindex(ep::ElectricPotential{T, N, S}, i::Int) where {T, N, S} = getindex(ep.data, i)
@inline getindex(ep::ElectricPotential{T, N, S}, s::Symbol) where {T, N, S} = getindex(ep.grid, s)


function ElectricPotential(fss::PotentialSimulationSetup{T, N, S} ; kwargs...)::ElectricPotential{T, N, S} where {T, N, S}
    return get_2π_potential(ElectricPotential{T, N, S}(fss.potential, fss.grid); kwargs...)
end


@recipe function f( ep::ElectricPotential{T, 3, :Cylindrical};
                    r = missing,
                    θ = missing,
                    z = missing,
                    contours_equal_potential=false ) where {T}
    g::Grid{T, 3, :Cylindrical} = ep.grid
   
    seriescolor --> :viridis
    st --> :heatmap
    aspect_ratio --> 1
    foreground_color_border --> nothing
    tick_direction --> :out
       
    cross_section::Symbol, idx::Int = if ismissing(θ) && ismissing(r) && ismissing(z)
        :θ, 1
    elseif !ismissing(θ) && ismissing(r) && ismissing(z)
        θ_rad::T = T(deg2rad(θ))
        while !(g[:θ].interval.left <= θ_rad <= g[:θ].interval.right)
            if θ_rad > g[:θ].interval.right
                θ_rad -= g[:θ].interval.right - g[:θ].interval.left
            elseif θ_rad < g[:θ].interval.left
                θ_rad += g[:θ].interval.right - g[:θ].interval.left
            end
        end
        :θ, searchsortednearest(g[:θ], θ_rad)
    elseif ismissing(θ) && !ismissing(r) && ismissing(z)
        :r, searchsortednearest(g[:r], T(r))
    elseif ismissing(θ) && ismissing(r) && !ismissing(z)
        :z, searchsortednearest(g[:z], T(z))
    else
        error(ArgumentError, ": Only one of the keywords `r, θ, z` is allowed.")
    end
    value::T = if cross_section == :θ
        g[:θ][idx]
    elseif cross_section == :r    
        g[:r][idx]
    elseif cross_section == :z
        g[:z][idx]
    end

    @series begin
        if cross_section == :θ
            title --> "Electric Potential @$(cross_section) = $(round(rad2deg(value), sigdigits = 2))"
            xlabel --> "r / m"
            ylabel --> "z / m"
            size --> ( 400, 350 / (g[:r][end] - g[:r][1]) * (g[:z][end] - g[:z][1]) )
            g[:r], g[:z], ep.data[:, idx,:]'
        elseif cross_section == :r
            title --> "Electric Potential @$(cross_section) = $(round(value, sigdigits = 2))"
            g[:θ], g[:z], ep.data[idx,:,:]'
        elseif cross_section == :z
            title --> "Electric Potential @$(cross_section) = $(round(value, sigdigits = 2))"
            proj --> :polar
            g[:θ], g[:r], ep.data[:,:,idx]
        end
    end
    if contours_equal_potential
        @series begin
            seriescolor := :thermal
            st := :contours
            if cross_section == :θ
                g[:r], g[:z], ep.data[:, idx,:]'
            elseif cross_section == :r
                g[:θ], g[:z], ep.data[idx,:,:]'
            elseif cross_section == :z
                proj --> :polar
                g[:θ], g[:r], ep.data[:,:,idx]
            end
        end
    end
end


function ElectricPotential(nt::NamedTuple)
    grid = Grid(nt.grid)
    T = typeof(ustrip(nt.values[1]))
    S = get_coordinate_type(grid)
    N = get_number_of_dimensions(grid)
    ElectricPotential{T, N, S}( ustrip.(uconvert.(u"V", nt.values)), grid)
end

Base.convert(T::Type{ElectricPotential}, x::NamedTuple) = T(x)

function NamedTuple(ep::ElectricPotential{T, 3, :Cylindrical}) where {T}
    return (
        grid = NamedTuple(ep.grid),
        values = ep.data * u"V",
    )
end

Base.convert(T::Type{NamedTuple}, x::ElectricPotential) = T(x)

