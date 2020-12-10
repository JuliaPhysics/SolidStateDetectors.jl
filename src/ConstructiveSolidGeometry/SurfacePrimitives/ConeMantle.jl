struct ConeMantle{T,TR,TP,TZ} <: AbstractSurfacePrimitive{T}
    r::TR #if not a Tuple, then ConeMantle is a Tube
    φ::TP
    z::TZ
    #if r is a Tuple, the first entry refers to the r at the bottom, the second one to the r at the top
    function ConeMantle( ::Type{T},
                   r::Union{T, Tuple{T,T}},
                   φ::Union{Nothing, <:AbstractInterval{T}},
                   z::Union{T, <:AbstractInterval{T}}) where {T, I<:AbstractInterval{T}}
        new{T,typeof(r),typeof(φ),typeof(z)}(r, φ, z)
    end
end

function ConeMantle(;rbot = 1, rtop = 0, φMin = 0, φMax = 2π, zMin = -1/2, zMax = 1/2)
    T = float(promote_type(typeof.((rbot, rtop, φMin, φMax, zMin, zMax))...))
    r = rbot == rtop ? T(rbot) : (T(rbot), T(rtop))
    φ = mod(T(φMax) - T(φMin), T(2π)) == 0 ? nothing : T(φMin)..T(φMax)
    z = zMax == -zMin ? T(zMax) : T(zMin)..T(zMax)
    ConeMantle( T, r, φ, z)
end
ConeMantle(rbot, rtop, φMin, φMax, zMin, zMax) = ConeMantle(;rbot, rtop, φMin, φMax, zMin, zMax)

function ConeMantle(rbot::R1, rtop::R2, height::H) where {R1<:Real, R2<:Real, H<:Real}
    T = float(promote_type(R1, R2, H))
    ConeMantle( T, (T(rbot), T(rtop)), nothing, T(height)/2)
end

get_r_limits(c::ConeMantle{T, T, <:Any, <:Any}) where {T} = (T(c.r), T(c.r))
get_r_limits(c::ConeMantle{T, <:Tuple, <:Any, <:Any}) where {T} = c.r

get_φ_limits(c::ConeMantle{T, <:Any, Nothing, <:Any}) where {T} = (T(0), T(2π), true)
get_φ_limits(c::ConeMantle{T, <:Any, <:AbstractInterval, <:Any}) where {T} = (c.φ.left, c.φ.right, false)

get_z_limits(c::ConeMantle{T}) where {T} = (_left_linear_interval(c.z), _right_linear_interval(c.z))

# plotting
function get_plot_points(c::ConeMantle{T}; n = 30) where {T <: AbstractFloat}

    plot_points = Vector{CartesianPoint{T}}[]

    rbot::T, rtop::T = get_r_limits(c)
    φMin::T, φMax::T, φ_is_full_2π::Bool = get_φ_limits(c)
    zMin::T, zMax::T = get_z_limits(c)

    φrange = range(φMin, φMax, length = n)

    #top and bottom circles
    rbot != 0 ? push!(plot_points, Vector{CartesianPoint{T}}([CartesianPoint{T}(rbot * cos(φ), rbot * sin(φ), zMin) for φ in φrange])) : nothing
    rtop != 0 ? push!(plot_points, Vector{CartesianPoint{T}}([CartesianPoint{T}(rtop * cos(φ), rtop * sin(φ), zMax) for φ in φrange])) : nothing

    #side line(s)
    for φ in (φ_is_full_2π ? T(0) : [φMin, φMax])
        if rbot != 0 || rtop != 0
            push!(plot_points, Vector{CartesianPoint{T}}([CartesianPoint{T}(rbot * cos(φ), rbot * sin(φ), zMin), CartesianPoint{T}(rtop * cos(φ), rtop * sin(φ), zMax)]))
        end
    end
    plot_points
end

function mesh(c::ConeMantle{T}; n = 30) where {T <: AbstractFloat}

    rbot::T, rtop::T = get_r_limits(c)
    φMin::T, φMax::T, φ_is_full_2π::Bool = get_φ_limits(c)
    zMin::T, zMax::T = get_z_limits(c)

    m = (rtop-rbot)/(zMax-zMin)
    φ = range(φMin, φMax, length = n+1)
    z = range(zMin, zMax, length = 2)

    X::Array{T,2} = [(m*(z[j]-zMin)+rbot)*cos(φ_i) for φ_i in φ, j in 1:length(z)]
    Y::Array{T,2} = [(m*(z[j]-zMin)+rbot)*sin(φ_i) for φ_i in φ, j in 1:length(z)]
    Z::Array{T,2} = [j for i in 1:length(φ), j in z]

    Mesh(X, Y, Z)
end
