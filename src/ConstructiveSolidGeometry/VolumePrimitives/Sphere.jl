struct Sphere{T,TR} <: AbstractVolumePrimitive{T}
    r::TR

    function Sphere( ::Type{T}, r::Union{T, <:AbstractInterval{T}}) where {T}
        new{T, typeof(r)}(r)
    end
end

#Constructors
function Sphere(; rMin = 0, rMax = 1)
    T = float(promote_type(typeof.((rMin, rMax))...))
    r = rMin == 0 ? T(rMax) : T(rMin)..T(rMax)
    Sphere(T, r)
end

Sphere(rMin, rMax) = Sphere(; rMin = rMin, rMax = rMax)
function Sphere(r::R) where {R <: Real}
    T = float(R)
    Sphere(T, T(r))
end

in(p::AbstractCoordinatePoint, s::Sphere) = _in_sph_r(p, s.r)

# read-in
function Geometry(::Type{T}, ::Type{Sphere}, dict::AbstractDict, input_units::NamedTuple) where {T}
    length_unit = input_units.length
    r = parse_r_of_primitive(T, dict, length_unit)
    return Sphere(T, r)
end

function Dictionary(s::Sphere{T}) where {T}
    dict = OrderedDict{String,Any}()
    dict["r"] = typeof(s.r) == T ? s.r : OrderedDict{String,Any}("from" => s.r.left, "to" => s.r.right)
    OrderedDict{String,Any}("sphere" => dict)
end

get_r_limits(s::Sphere{T}) where {T} = (_left_radial_interval(s.r), _right_radial_interval(s.r))

get_decomposed_surfaces(s::Sphere{T, T}) where {T} = AbstractSurfacePrimitive[SphereMantle{T}(s.r)]

function get_decomposed_surfaces(s::Sphere{T, <:AbstractInterval{T}}) where {T}
    rMin::T, rMax::T = get_r_limits(s)
    AbstractSurfacePrimitive[SphereMantle{T}(rMin), SphereMantle{T}(rMax)]
end

function sample(s::Sphere{T}, Nsamps::NTuple{3,Int} = (2,5,3))::Vector{CylindricalPoint{T}} where {T}
    rMin::T, rMax::T = get_r_limits(s)
    samples = [
        CylindricalPoint{T}(r,φ,z)
        for z in (Nsamps[3] ≤ 1 ? rMax : range(-rMax, rMax, length = Nsamps[3]))
        for r in (Nsamps[1] ≤ 1 ? rMax : range((abs(z) > rMin ? 0 : sqrt(rMin^2 - z^2)), sqrt(rMax^2 - z^2), length = Nsamps[1]))
        for φ in (Nsamps[2] ≤ 1 ? 0 : range(0, 2π, length = Nsamps[2]))
    ]
end
