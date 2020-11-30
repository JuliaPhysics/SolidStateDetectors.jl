struct ScaledGeometry{T,P<:AbstractGeometry{T}} <: AbstractGeometry{T}
    p::P
    inv_s::SVector{3,T}
    ScaledGeometry(p::AbstractGeometry{T}, s::SVector{3,T}) where {T} = new{T,typeof(p)}(p, inv.(s))
    ScaledGeometry(g::ScaledGeometry{T}, s::SVector{3,T}) where {T} = new{T,typeof(g.p)}(g.p, g.inv_s .* inv.(s))
end
in(p::CartesianPoint, g::ScaledGeometry) = in(CartesianPoint(g.inv_s .* p), g.p)
scale(g::AbstractGeometry, s::SVector{3,T}) where {T} = ScaledGeometry(g, s)

struct RotatedGeometry{T,P<:AbstractGeometry{T},RT} <: AbstractGeometry{T}
    p::P
    inv_r::RotMatrix{3,RT,9}
    RotatedGeometry(p::AbstractGeometry{T}, r::RotMatrix{3,RT,9}) where {T,RT} = new{T,typeof(p),RT}(p, inv(r))
    RotatedGeometry(g::RotatedGeometry{T}, r::RotMatrix{3,RT,9}) where {T,RT} = new{T,typeof(g.p),RT}(g.p, g.inv_r * inv(r))
end
in(p::CartesianPoint, g::RotatedGeometry) = in(g.inv_r * p, g.p)
rotate(g::AbstractGeometry{T}, r::RotMatrix{3,RT,9}) where {T,RT} = RotatedGeometry(g, r)

struct TranslatedGeometry{T,P<:AbstractGeometry{T}} <: AbstractGeometry{T}
    p::P
    t::CartesianVector{T}
end
in(p::CartesianPoint, g::TranslatedGeometry) = in(p - g.t, g.p)
translate(g::AbstractGeometry{T}, t::CartesianVector{T}) where {T} = TranslatedGeometry(g, t)
translate(g::TranslatedGeometry{T}, t::CartesianVector{T}) where {T} = TranslatedGeometry(g.p, g.t + t)

