struct ScaledGeometry{T,P<:AbstractGeometry{T}} <: AbstractGeometry{T}
    p::P
    inv_s::SVector{3,T}
    ScaledGeometry(p::P, s::SVector{3,T}) where {T,P} = new{T,P}(p, inv.(s))
end
in(p::CartesianPoint, g::ScaledGeometry) = in(CartesianPoint(g.inv_s .* p), g.p)
scale(g::AbstractGeometry, s::SVector{3,T}) where {T} = ScaledGeometry(g, s)
scale(g::ScaledGeometry, s::SVector{3,T}) where {T} = ScaledGeometry(g.p, inv.(g.inv_s) .* s)

struct RotatedGeometry{T,P<:AbstractGeometry{T},RT} <: AbstractGeometry{T}
    p::P
    inv_r::RotMatrix{3,RT,9}
    RotatedGeometry(p::AbstractGeometry{T}, r::RotMatrix{3,RT,9}) where {T,RT} = new{T,typeof(p),RT}(p, inv(r))
end
in(p::CartesianPoint, g::RotatedGeometry) = in(g.inv_r * p, g.p)
rotate(g::AbstractGeometry{T}, r::RotMatrix{3,RT,9}) where {T,RT} = RotatedGeometry(g, r)
rotate(g::RotatedGeometry{T,<:Any,RT}, r::RotMatrix{3,RT,9}) where {T,RT} = RotatedGeometry(g.p, r * inv(g.inv_r))


struct TranslatedGeometry{T,P<:AbstractGeometry{T}} <: AbstractGeometry{T}
    p::P
    t::CartesianVector{T}
end
in(p::CartesianPoint, g::TranslatedGeometry) = in(p - g.t, g.p)
translate(g::AbstractGeometry{T}, t::CartesianVector{T}) where {T} = TranslatedGeometry(g, t)
translate(g::TranslatedGeometry{T}, t::CartesianVector{T}) where {T} = TranslatedGeometry(g.p, g.t + t)

