csg_round_lin(x::T; sigdigits::Int = 9) where {T} = round(x, sigdigits = sigdigits) # min. Δx = 1 nm 
csg_round_lin(x::Float32; sigdigits::Int = 5) where {T} = round(x, sigdigits = sigdigits) # min. Δx = 10 μm 
csg_round_rad(φ::T; digits::Int = 4) where {T} = round(φ, digits = digits) # min. Δφ ≈ 0.011 °

csg_default_tol(::Type{T}) where {T} = T(1e-12)

csg_isapprox(x::T, y::T; csgtol::T = T(1e-9)) where {T} = abs(x - y) < csgtol

