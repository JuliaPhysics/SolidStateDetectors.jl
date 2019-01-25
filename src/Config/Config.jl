"""
    abstract type AbstractConfig{T <: AbstractFloat} end

`T`: Type of points or precision. 
"""
abstract type AbstractConfig{T <: AbstractFloat} end

function get_ρ_and_ϵ(pt, config::AbstractConfig{T})::Tuple{T, T} where {T <: AbstractFloat}
    error("User has to define a method `get_ρ_and_ϵ` for `pt::$(typeof(pt))` and `config::$(typeof(config))`.")
end
function set_pointtypes_and_fixed_potentials!(pointtypes::Array{PointType, N}, potential::Array{T, N}, grid::Grid{T, N, S}, config::AbstractConfig{T})::Nothing where {T, N, G, S} 
    error("User has to define a method `set_pointtypes_and_fixed_potentials!` for `config::$(typeof(config))`.")
end



function get_ρ_and_ϵ(pt::Cylindrical{T}, ssd::SolidStateDetector{T})::Tuple{T, T} where {T <: AbstractFloat}
    if in(pt, ssd)
        ρ::T = get_charge_density(ssd, pt.r, pt.θ, pt.z) * elementary_charge
        ϵ::T = ssd.material_detector.ϵ_r
        return ρ, ϵ
    else
        ρ = 0
        ϵ = ssd.material_environment.ϵ_r 
        return ρ, ϵ
    end
end
function set_pointtypes_and_fixed_potentials!(pointtypes::Array{PointType, N}, potential::Array{T, N}, 
        grid::Grid{T, N, :Cylindrical}, ssd::SolidStateDetector{T}; weighting_potential_channel_idx::Union{Missing, Int} = missing)::Nothing where {T <: AbstractFloat, N}
    
    channels::Array{Int, 1} = if !ismissing(weighting_potential_channel_idx)
        ssd.grouped_channels[weighting_potential_channel_idx]
    else
        Int[]
    end

    axr::Vector{T} = grid[:r].ticks
    axθ::Vector{T} = grid[:θ].ticks
    axz::Vector{T} = grid[:z].ticks
    for iz in axes(potential, 3)
        z::T = axz[iz]
        for iθ in axes(potential, 2)
            θ::T = axθ[iθ]
            for ir in axes(potential, 1)
                r::T = axr[ir]
                pt::Cylindrical{T} = Cylindrical{T}( r, θ, z )              

                if is_boundary_point(ssd, r, θ, z, axr, axθ, axz)
                    pot::T = if ismissing(weighting_potential_channel_idx)
                        get_boundary_value( ssd, r, θ, z, axr)
                    else
                        in(get_segment_idx(ssd, r, θ, z), channels) ? 1 : 0
                    end
                    potential[ ir, iθ, iz ] = pot
                    pointtypes[ ir, iθ, iz ] = zero(PointType)
                elseif in(pt, ssd)
                    pointtypes[ ir, iθ, iz ] += pn_junction_bit 
                end

            end
        end
    end
    nothing
end

