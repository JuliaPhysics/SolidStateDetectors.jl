struct Sphere{T,TR} <: AbstractVolumePrimitive{T}
    r::TR
    
    function Sphere( ::Type{T}, r::Union{T, <: AbstractInterval{T}}) where {T}
        new{T, typeof(r)}(r)
    end
end

#Constructors
function Sphere(; rMin = 0, rMax = 1)
    T = float(promote_type(typeof.((rMin, rMax))...))
    r = rMin == 0 ? T(rMax) : T(rMin)..T(rMax)
    Sphere(T, r)
end

Sphere(rMin, rMax) = Sphere(;rMin, rMax)
function Sphere(r::R) where {R <: Real}
    T = float(R)
    Sphere(T, T(r))
end

in(p::AbstractCoordinatePoint, s::Sphere{<:Any, <:Any}) = _in_sph_r(p, s.r)


# read-in
function Geometry(T::DataType, t::Val{:sphere}, dict::Union{Dict{String,Any}, Dict{Any,Any}}, iud::Dict{String,Unitful.Units})
    length_unit = iud["length"]
    r = _get_r_of_primitive(T, dict["r"], length_unit)
    return Sphere(T, r)
end


# plotting
function get_plot_points(s::Sphere{T}) where {T <: AbstractFloat}
    
    plot_points = Vector{CartesianPoint{T}}[]
    
    rMin::T = _left_radial_interval(s.r)
    rMax::T = _right_radial_interval(s.r)
    φrange = range(0, 2π, length = 36)
    
    for r in (rMin == 0 ? [rMax] : [rMin, rMax])
        for φ in range(0, 2π, length = 11)
        push!(plot_points, Vector{CartesianPoint{T}}([CartesianPoint{T}(r * sin(θ) * cos(φ), r * sin(θ) * sin(φ), r * cos(θ)) for θ in φrange]))
        end
        push!(plot_points, Vector{CartesianPoint{T}}([CartesianPoint{T}(r * cos(φ), r * sin(φ), 0) for φ in φrange]))
    end

    plot_points
end