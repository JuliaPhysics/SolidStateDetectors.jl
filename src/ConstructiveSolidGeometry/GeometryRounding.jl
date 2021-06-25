csg_round_lin(x::T; sigdigits::Int = 9) where {T} = round(x, sigdigits = sigdigits) # min. Δx = 1 nm 
csg_round_lin(x::Float32; sigdigits::Int = 5) where {T} = round(x, sigdigits = sigdigits) # min. Δx = 10 μm 
csg_round_rad(φ::T; digits::Int = 4) where {T} = round(φ, digits = digits) # min. Δφ ≈ 0.011 °

csg_isapprox_lin(x::T, y::T; Δ::T = T(1e-9)) where {T} = abs(x - y) < Δ
csg_isapprox_rad(α::T, β::T; Δ::T = T(1e-4)) where {T} = abs(α - β) < Δ