@inline _linear_endpoints(z::Real) = (-z, z)
@inline _linear_endpoints(z::AbstractInterval) = endpoints(z)
@inline _linear_endpoints(r::Tuple{T,T}) where {T} = r

@inline _radial_endpoints(r::T) where {T <: Real} = (zero(T), r)
@inline _radial_endpoints(r::AbstractInterval) = endpoints(r)
@inline _radial_endpoints(r::Tuple{T,T}) where {T} = r

@inline _in_angular_interval_closed(α::T, α_int::T; csgtol::T = csg_default_tol(T)) where {T} = mod(α + csgtol, T(2π)) ≤ α_int + csgtol
@inline function _in_angular_interval_closed(α::T, α_int::Tuple{T,T}; csgtol::T = csg_default_tol(T)) where {T} 
    m = mod(α - α_int[1], T(2π))
    d = (α_int[2] - α_int[1])
    m ≤ d + csgtol 
end

@inline _in_angular_interval_open(α::T, α_int::T; csgtol::T = csg_default_tol(T)) where {T} = csgtol < mod(α, T(2π)) ≤ α_int - csgtol
@inline function _in_angular_interval_open(α::T, α_int::Tuple{T,T}; csgtol::T = csg_default_tol(T)) where {T} 
    csgtol < mod(α - α_int[1], T(2π)) < (α_int[2] - α_int[1]) - csgtol 
end

@inline _in_angular_interval_closed(α::Real, α_int::Nothing, tol::Real = 0) = true
# @inline _in_angular_interval_closed(α::Real, α_int::AbstractInterval{T}, tol::Real = 0) where {T} = mod(α - (α_int.left-tol), T(2π)) ≤ (α_int.right+tol) - (α_int.left-tol)
# 
# @inline _in_angular_interval_open(α::T, α_int::Nothing) where {T<:Real} = 0 < mod(α, T(2π)) < T(2π)
# @inline _in_angular_interval_open(α::Real, α_int::AbstractInterval{T}) where {T} = 0 < mod(α - α_int.left, T(2π)) < width(α_int)

@inline _in_φ(p::CylindricalPoint, φ::Nothing) = true
@inline function _in_φ(p::CylindricalPoint, φ::Real)
    φp = mod(float(p.φ), 2π)
    φmax = mod(float(φ), 2π)
    return φp <= φmax
end
@inline _in_φ(p::CartesianPoint, φ::Union{Nothing, Real}) =  _in_φ(CylindricalPoint(p), φ)

"""
    _φNear(φ, φMax)
Return the nearest azimuthal boundary (0 or φMax) for segmented surfaces with φ ∈ [0, φMax].
"""
@inline function _φNear(φ::T, φMax::T) where T
    φ = mod(φ, 2π)
    φMax = mod(φMax, 2π)
    # distance to lower boundary (0)
    d0 = min(φ, 2π - φ)
    # distance to upper boundary (φMax)
    dMax = min(abs(φ - φMax), 2π - abs(φ - φMax))
    
    return d0 <= dMax ? zero(T) : φMax
end
