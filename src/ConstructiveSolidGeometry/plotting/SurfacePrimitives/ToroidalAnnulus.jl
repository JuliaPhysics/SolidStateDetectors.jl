function get_plot_points(t::ToroidalAnnulus{T}; n = 30) where {T <: AbstractFloat}

    plot_points = Vector{CartesianPoint{T}}[]

    r_tubeMin::T, r_tubeMax::T = get_r_tube_limits(t)
    θMin::T, θMax::T, θ_is_full_2π::Bool = get_θ_limits(t)
    θrange = range(θMin, θMax, length = n + 1)

    #circle(s)
    for r_tube in [r_tubeMin, r_tubeMax]
        if r_tube == 0 continue end
        push!(plot_points, Vector{CartesianPoint{T}}([CartesianPoint{T}((t.r_torus+r_tube*cos(θ))cos(t.φ), (t.r_torus+r_tube*cos(θ))sin(t.φ), r_tube*sin(θ)) for θ in θrange]))
    end

    #for incomplete θ: lines of cross-sections
    if !θ_is_full_2π
        for θ in [θMin, θMax]
            push!(plot_points, Vector{CartesianPoint{T}}([CartesianPoint{T}((t.r_torus+r_tubeMin*cos(θ))cos(t.φ), (t.r_torus+r_tubeMin*cos(θ))sin(t.φ), r_tubeMin*sin(θ)), CartesianPoint{T}((t.r_torus+r_tubeMax*cos(θ))cos(t.φ), (t.r_torus+r_tubeMax*cos(θ))sin(t.φ), r_tubeMax*sin(θ))]))
        end
    end
    plot_points
end

function mesh(t::ToroidalAnnulus{T}; n = 30) where {T <: AbstractFloat}

    r_tubeMin::T, r_tubeMax::T = get_r_tube_limits(t)
    θMin::T, θMax::T, _ = get_θ_limits(t)
    r_tube = range(r_tubeMin, r_tubeMax, length = 2)
    sφ, cφ = sincos(t.φ)
    θrange = range(θMin, θMax, length = n + 1)
    sθrange = sin.(θrange)
    cθrange = cos.(θrange)

    X::Array{T,2} = [(t.r_torus+r*cθ)*cφ for cθ in cθrange, r in r_tube]
    Y::Array{T,2} = [(t.r_torus+r*cθ)*sφ for cθ in cθrange, r in r_tube]
    Z::Array{T,2} = [r*sθ for sθ in sθrange, r in r_tube]

    Mesh(X, Y, Z)
end
