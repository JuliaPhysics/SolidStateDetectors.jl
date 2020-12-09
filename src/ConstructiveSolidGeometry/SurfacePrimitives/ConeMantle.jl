struct ConeMantle{T,TR,TP,TZ} <: AbstractSurfacePrimitive{T}
    r::TR #if not a Tuple, then Cone is a Tube
    φ::TP
    z::TZ
    #if r is a Tuple, the first entry refers to the r-interval at the bottom, the second one to the r-interval at the top
    function ConeMantle( ::Type{T},
                   r::Union{T, Tuple{T,T}},
                   φ::Union{Nothing, <:AbstractInterval{T}},
                   z::Union{T, <:AbstractInterval{T}}) where {T, I<:AbstractInterval{T}}
        new{T,typeof(r),typeof(φ),typeof(z)}(r, φ, z)
    end
end

function ConeMantle(;rbot = 0, rtop = 0, φMin = 0, φMax = 2π, zMin = -1/2, zMax = 1/2)
    T = float(promote_type(typeof.((rtop, rbot, φMin, φMax, zMin, zMax))...))
    r = rbot == rtop ? rbot : (T(rbot), T(rtop))
    φ = mod(T(φMax) - T(φMin), T(2π)) == 0 ? nothing : T(φMin)..T(φMax)
    z = zMax == -zMin ? T(zMax) : T(zMin)..T(zMax)
    Cone( T, r, φ, z)
end
ConeMantle(rbotMin, rbotMax, rtopMin, rtopMax, φMin, φMax, zMin, zMax) = Cone(;rbotMin, rbotMax, rtopMin, rtopMax, φMin, φMax, zMin, zMax)

function ConeMantle(rbot::R1, rtop::R2, height::H) where {R1<:Real, R2<:Real, H<:Real}
    T = float(promote_type(R1, R2, H))
    Cone( T, (T(rbot), T(rtop)), nothing, T(height)/2)
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
