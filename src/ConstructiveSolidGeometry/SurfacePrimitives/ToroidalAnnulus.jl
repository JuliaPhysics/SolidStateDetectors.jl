struct ToroidalAnnulus{T,TR,TB,TP,TT} <: AbstractSurfacePrimitive{T}
    r_torus::TR
    r_tube::TB
    φ::TP
    θ::TT
    function ToroidalAnnulus( ::Type{T},
                   r_torus::T,
                   r_tube::Union{T, <:AbstractInterval{T}},
                   φ::T,
                   θ::Union{Nothing, <:AbstractInterval{T}}) where {T}
        new{T,T,typeof(r_tube),T,typeof(θ)}(r_torus, r_tube, φ, θ)
    end
end

#Constructors
function ToroidalAnnulus(;r_torus = 1, r_tubeMin = 0, r_tubeMax = 1, φ = 0, θMin = 0, θMax = 2π)
    T = float(promote_type(typeof.((r_torus, r_tubeMin, r_tubeMax, φ, θMin, θMax))...))
    r_tube = r_tubeMin == 0 ? T(r_tubeMax) : T(r_tubeMin)..T(r_tubeMax)
    θ = mod(T(θMax) - T(θMin), T(2π)) == 0 ? nothing : T(θMin)..T(θMax)
    ToroidalAnnulus( T, T(r_torus), r_tube, T(φ), θ)
end
ToroidalAnnulus(r_torus, r_tubeMin, r_tubeMax, φ, θMin, θMax) = ToroidalAnnulus(;r_torus = r_torus, r_tubeMin = r_tubeMin, r_tubeMax = r_tubeMax, φ = φ, θMin = θMin, θMax = θMax)

#use _in_torr_r_tube functions from Torus, these will have to be moved to PointsAndVectors
