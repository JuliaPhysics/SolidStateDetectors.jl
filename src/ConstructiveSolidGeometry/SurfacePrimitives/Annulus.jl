struct Annulus{T,TR,TP} <: AbstractSurfacePrimitive{T}
    r::TR
    φ::TP
    function Annulus( ::Type{T},
                   r::Union{T, <:AbstractInterval{T}},
                   φ::Union{Nothing, <:AbstractInterval{T}}) where {T}
        new{T,typeof(r),typeof(φ)}(r, φ)
    end
end

#Constructors
function Annulus(; rMin = 0, rMax = 1, φMin = 0, φMax = 2π)
    T = float(promote_type(typeof.((rMin, rMax, φMin, φMax))...))
    r = rMin == 0 ? T(rMax) : T(rMin)..T(rMax)
    φ = mod(T(φMax) - T(φMin), T(2π)) == 0 ? nothing : T(φMin)..T(φMax)
    Annulus(T, r, φ)
end

Annulus(rMin, rMax, φMin, φMax) = Annulus(;rMin, rMax, φMin, φMax)

function Annulus(r::R) where {R <: Real}
    T = float(R)
    Annulus(T, T(r), nothing)
end

get_r_limits(a::Annulus{T, <:Union{T, AbstractInterval{T}}, <:Any}) where {T} =
    (_left_radial_interval(a.r), _right_radial_interval(a.r))

get_φ_limits(a::Annulus{T, <:Any, Nothing}) where {T} = (T(0), T(2π), true)
get_φ_limits(a::Annulus{T, <:Any, <:AbstractInterval}) where {T} = (a.φ.left, a.φ.right, false)

get_z_limits(a::Annulus{T}) where {T} = nothing
get_θ_limits(a::Annulus{T}) where {T} = nothing

#plotting
function get_plot_points(a::Annulus{T}; n = 30) where {T <: AbstractFloat}

    plot_points = Vector{CartesianPoint{T}}[]

    rMin::T, rMax::T = get_r_limits(a)
    φMin::T, φMax::T, φ_is_full_2π::Bool = get_φ_limits(a)
    φrange = range(φMin, φMax, length = n)
    z::T = 0

    #circle(s)
    for r in [rMin, rMax]
        if r == 0 continue end
        push!(plot_points, Vector{CartesianPoint{T}}([CartesianPoint{T}(r * cos(φ), r * sin(φ), z) for φ in φrange]))
    end

    #for incomplete φ: lines of cross-sections
    if !φ_is_full_2π
        for φ in [φMin, φMax]
            push!(plot_points, Vector{CartesianPoint{T}}([CartesianPoint{T}(rMin * cos(φ), rMin * sin(φ), z), CartesianPoint{T}(rMax * cos(φ), rMax * sin(φ), z)]))
        end
    end
    plot_points
end

function mesh(a::Annulus{T}; n = 30) where {T <: AbstractFloat}
    rMin::T, rMax::T = get_r_limits(a)
    φMin::T, φMax::T, φ_is_full_2π::Bool = get_φ_limits(a)
    φ = range(φMin, φMax, length = n+1)
    r = range(rMin, rMax, length = 2)
    z = zeros(length(r))
    X::Array{T,2} = [r[j]*cos(φ_i) for φ_i in φ, j in 1:length(r)]
    Y::Array{T,2} = [r[j]*sin(φ_i) for φ_i in φ, j in 1:length(r)]
    Z::Array{T,2} = [j for i in 1:length(φ), j in z]
    Mesh(X, Y, Z)
end
