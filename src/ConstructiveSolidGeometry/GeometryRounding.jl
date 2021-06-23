csg_round_lin(x::T) where {T} = round(x, sigdigits = 9) # min. Δx = 1 nm 
csg_round_lin(x::Float32) where {T} = round(x, sigdigits = 5) # min. Δx = 10 μm 
csg_round_rad(φ::T) where {T} = round(φ, digits = 4) # min. Δφ ≈ 0.011 °

csg_isapprox_lin(x::T, y::T) where {T} = csg_round_lin(x) == csg_round_lin(y) #abs(x - y) < T(1e-9)
csg_isapprox_rad(α::T, β::T) where {T} = csg_round_rad(α) == csg_round_rad(β) #abs(x - y) < T(1e-4)