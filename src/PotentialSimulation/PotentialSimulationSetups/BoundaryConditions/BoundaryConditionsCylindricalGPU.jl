

function apply_boundary_conditions!(pssrb::PotentialSimulationSetupRBGPU{T, 3, 4, Cylindrical, AT}, update_even_points, only2d) where {T, AT}
    rbi::Int = update_even_points ? rb_even::Int : rb_odd::Int
    # nrbi::Int = update_even_points ? rb_odd::Int : rb_even::Int
    p = pssrb.potential
    @inbounds p[:, :, end, rbi] .= pssrb.grid_boundary_factors[1][2] .* view(p, :, :, size(p, 3) - 2, rbi) # infinity boundaries in r        
    if !only2d
        # # # Azimutal-Axis
        @inbounds p[:,   1, :, rbi] .= view(p, :,       3, :, rbi) # reflecting in phi
        @inbounds p[:, end, :, rbi] .= view(p, :, size(p, 2) - 2, :, rbi) # reflecting in phi
    end
    @inbounds p[  1, :, :, rbi] .= pssrb.grid_boundary_factors[3][1] .* view(p,       2, :, :, rbi)  # infinity boundaries in z (only ±1 because this is the compressed (redblack) dimension)
    @inbounds p[end, :, :, rbi] .= pssrb.grid_boundary_factors[3][2] .* view(p, size(p, 1) - 1, :, :, rbi)    # infinity boundaries in z (only ±1 because this is the compressed (redblack) dimension)
    nothing
end