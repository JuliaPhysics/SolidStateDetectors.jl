@inline function sor_kernel(
    potential::AbstractArray{T, 4},
    imp_scale::AbstractArray{T, 4},
    point_types::AbstractArray{PointType, 4},
    volume_weights::AbstractArray{T, 4},
    q_eff_imp::AbstractArray{T, 4},
    q_eff_fix::AbstractArray{T, 4},
    ϵ_r::AbstractArray{T, 3},
    geom_weights::NTuple{3, AbstractArray{T, 2}},
    sor_const::AbstractArray{T, 1},
    update_even_points::Bool,
    depletion_handling_enabled::Bool,
    is_weighting_potential::Bool,
    only2d::Bool,
    ::Type{S},
    gpu_inds::NTuple{3, Int}
) where {T, S}
    # Comparison to CPU indices: (Cyl / Car)
    # i3 <-> idx3 / ir / iz
    # i2 <-> idx2 / iφ / iy
    # i1 <-> idx1 / iz / ix
    i1, in2, in3 = gpu_inds
    i1 += 1
    i2 = in2 + 1
    i3 = in3 + 1
    in1 = nidx(i1, update_even_points, iseven(i2 + i3))
    i1r = get_rbidx_right_neighbour(i1, update_even_points, iseven(i2 + i3))
    
    rb_tar_idx, rb_src_idx = update_even_points ? (rb_even::Int, rb_odd::Int) : (rb_odd::Int, rb_even::Int) 

    geom_weights_3 = get_geom_weights_outerloop(geom_weights, in3, S)
    geom_weights_2 = prepare_weights_in_middleloop(
        geom_weights, S, i2, in2, 
        geom_weights_3...,
        in3 == 1
    )
    weights = calculate_sor_weights(
        in1, 
        S, 
        ϵ_r,
        geom_weights[3],
        i2, in2, i3, in3,
        geom_weights_2...
    )

    old_potential = potential[i1, i2, i3, rb_tar_idx]
    q_eff = is_weighting_potential ? zero(T) : (q_eff_imp[i1, i2, i3, rb_tar_idx] * imp_scale[i1, i2, i3, rb_tar_idx] + q_eff_fix[i1, i2, i3, rb_tar_idx])

    neighbor_potentials = get_neighbor_potentials(
        potential, old_potential, i1, i2, i3, i1r, in2, in3, rb_src_idx, only2d
    )      
    
    new_potential = calc_new_potential_SOR_3D(
        q_eff,
        volume_weights[i1, i2, i3, rb_tar_idx],
        weights,
        neighbor_potentials,
        old_potential,
        get_sor_constant(sor_const, S, in3)
    )

    if depletion_handling_enabled
        new_potential, imp_scale[i1, i2, i3, rb_tar_idx] = handle_depletion(
            new_potential,
            imp_scale[i1, i2, i3, rb_tar_idx],
            r0_handling_depletion_handling(neighbor_potentials, S, in3),
            q_eff_imp[i1, i2, i3, rb_tar_idx],
            volume_weights[i1, i2, i3, rb_tar_idx],
            get_sor_constant(sor_const, S, in3)
        )
    end

    potential[i1, i2, i3, rb_tar_idx] = ifelse(point_types[i1, i2, i3, rb_tar_idx] & update_bit > 0, new_potential, old_potential)
    nothing
end

@inline function get_neighbor_potentials(
    potential::AbstractArray{T, 4},
    old_potential, i1, i2, i3, i1r, in2, in3, rb_src_idx, only2d::Bool
)::NTuple{6, T} where {T}
    @inbounds return ( # p: potential; 1: RB-dimension; l/r: left/right
    potential[i1r - 1,     i2,     i3, rb_src_idx], # p1l
    potential[    i1r,     i2,     i3, rb_src_idx], # p1r
    only2d ? old_potential : potential[ i1,    in2, i3, rb_src_idx], # p2l
    only2d ? old_potential : potential[ i1, i2 + 1, i3, rb_src_idx], # p2r
    potential[     i1,     i2,    in3, rb_src_idx], # p3l
    potential[     i1,     i2, i3 + 1, rb_src_idx], # p3r
    ) 
end

@inline function prepare_weights_in_middleloop(
    geom_weights::NTuple{3, <:AbstractArray{T, 2}}, ::Type{Cylindrical},
    i2, in2,
    pwwrr, pwwrl, r_inv_pwΔmpr, Δr_ext_inv_r_pwmprr, Δr_ext_inv_l_pwmprl, Δmpr_squared, 
    is_r0::Bool
) where {T}
    pwwφr        = geom_weights[2][1, in2]
    pwwφl        = geom_weights[2][2, in2]
    pwΔmpφ       = geom_weights[2][3, in2]
    Δφ_ext_inv_r = geom_weights[2][4,  i2]
    Δφ_ext_inv_l = geom_weights[2][4, in2]

    if is_r0
        pwwφr = T(0.5)
        pwwφl = T(0.5)
        pwΔmpφ = T(2π)
        Δφ_ext_inv_r = inv(pwΔmpφ)
        Δφ_ext_inv_l = Δφ_ext_inv_r
    end
    pwwrr_pwwφr = pwwrr * pwwφr
    pwwrl_pwwφr = pwwrl * pwwφr
    pwwrr_pwwφl = pwwrr * pwwφl
    pwwrl_pwwφl = pwwrl * pwwφl

    pwΔmpφ_Δmpr_squared = pwΔmpφ * Δmpr_squared
    Δr_ext_inv_r_pwmprr_pwΔmpφ = Δr_ext_inv_r_pwmprr * pwΔmpφ
    Δr_ext_inv_l_pwmprl_pwΔmpφ = Δr_ext_inv_l_pwmprl * pwΔmpφ
    r_inv_pwΔmpr_Δφ_ext_inv_r = r_inv_pwΔmpr * Δφ_ext_inv_r
    r_inv_pwΔmpr_Δφ_ext_inv_l = r_inv_pwΔmpr * Δφ_ext_inv_l
    return (
        pwwrr, pwwrl, pwwφr, pwwφl, 
        pwwrr_pwwφr, pwwrl_pwwφr, pwwrr_pwwφl, pwwrl_pwwφl,
        pwΔmpφ_Δmpr_squared,
        Δr_ext_inv_r_pwmprr_pwΔmpφ, Δr_ext_inv_l_pwmprl_pwΔmpφ,
        r_inv_pwΔmpr_Δφ_ext_inv_r, r_inv_pwΔmpr_Δφ_ext_inv_l
    )
end

"""
    function get_sor_kernel(::Type{S}, args...) 

where S is either Cartesian or Cylindrical. 

Developer notes:
Currently (February 2022), there are some limitations to the `@kernel` macro 
of the package KernelAbstractions.jl. Especially, regarding usage of dispatch. 

Thus, we have to write two kernel functions right now for the Cartesian & Cylindrical case:
`sor_cyl_gpu!` and `sor_car_gpu!`.

Inside kernel functions, everything is (and has to be) inlined and we can make use of multiple dispatch. 
So in the end we only have to write one function for the kernel, `sor_kernel`, 
which is then inlined inside the two kernel functions.

We can also use most of the CPU functions with the restriction that all 
types have to be independent on the GPU-indices of the kernel. 
E.g., making use of `i23_is_even_t = Val(iseven(i2 + i3))` and other similar statements is not possible,
which are used in the CPU implementation for optimization.
On the GPU (currently) those statements have to be calculated and booleans have to be passed.
Maybe this will change in the future.
"""
get_sor_kernel(::Type{Cylindrical}, args...) = sor_cyl_gpu!(args...)
get_sor_kernel(::Type{Cartesian},   args...) = sor_car_gpu!(args...)
get_sor_kernel(::Type{Cylindrical}, ::CPU) = nothing
get_sor_kernel(::Type{Cartesian}, ::CPU) = nothing

get_device(::Type{Array}) = CPU() 
