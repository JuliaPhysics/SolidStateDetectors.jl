struct StrechedGeometry{T,P} <: AbstractGeometry{T}
    p::P
    inv_s::SVector{3,T}
    StrechedGeometry(p::P, s::SVector{3,T}) where {T,P} = new{T,P}(p, inv.(s))
end
in(p::CartesianPoint, sp::StrechedGeometry) = in(CartesianPoint(sp.inv_s .* p), sp.p)

struct RotatedGeometry{T,P,RT} <: AbstractGeometry{T}
    p::P
    r::RotMatrix{3,RT,9}
end
in(p::CartesianPoint, rp::RotatedGeometry) = in(rp.r * p, rp.p)

struct TranslatedGeometry{T,P} <: AbstractGeometry{T}
    p::P
    t::CartesianVector{T}
end
in(p::CartesianPoint, tp::TranslatedGeometry) = in(p - tp.t, tp.p)



