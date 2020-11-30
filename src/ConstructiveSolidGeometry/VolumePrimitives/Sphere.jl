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