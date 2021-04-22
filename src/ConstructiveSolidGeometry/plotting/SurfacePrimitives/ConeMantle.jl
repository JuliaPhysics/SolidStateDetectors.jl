function get_plot_points(c::ConeMantle{T}; n = 30) where {T <: AbstractFloat}

    plot_points = Vector{CartesianPoint{T}}[]

    rbot::T, rtop::T = get_r_limits(c)
    φMin::T, φMax::T, φ_is_full_2π::Bool = get_φ_limits(c)
    zMin::T, zMax::T = get_z_limits(c)

    φrange = range(φMin, φMax, length = n + 1)

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
    f = (φMax - φMin)/(2π)
    n = Int(ceil(n*f))
    m = (rtop-rbot)/(zMax-zMin)
    φ = range(φMin, φMax, length = n+1)
    z = range(zMin, zMax, length = 2)

    X::Array{T,2} = [(m*(z[j]-zMin)+rbot)*cos(φ_i) for φ_i in φ, j in 1:length(z)]
    Y::Array{T,2} = [(m*(z[j]-zMin)+rbot)*sin(φ_i) for φ_i in φ, j in 1:length(z)]
    Z::Array{T,2} = [j for i in 1:length(φ), j in z]

    Mesh(X, Y, Z)
end
