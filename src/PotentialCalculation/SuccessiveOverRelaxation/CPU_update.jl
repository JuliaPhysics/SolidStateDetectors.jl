"""
    function update!(   
        pcs::PotentialCalculationSetup{T}; 
        ::Nothing, # these two unused arguments are used such that the method
        ::Any;     # is similar to the GPU method for it.
        use_nthreads::Int = Base.Threads.nthreads(), 
        depletion_handling::Val{depletion_handling_enabled} = Val{false}(), 
        is_weighting_potential::Val{_is_weighting_potential} = Val{false}(),
        only2d::Val{only_2d} = Val{false}()
)::Nothing where {T, depletion_handling_enabled, only_2d, _is_weighting_potential}

This function performs one iteration of the SOR. One iteration consists out of 4 steps:

    1) Iterate in parallel over all even points and update their potential. 
    2) Apply the boundary conditions at the ends of the grid for all even points. 
    3) Iterate in parallel over all odd points and update their potential. 
    2) Apply the boundary conditions at the ends of the grid for all odd points. 
"""
@inline function update!(   
    pcs::PotentialCalculationSetup{T},
    ::Nothing, # these to unused arguments are used such that the method
    ::Any;     # is similar to the GPU method for it.
    use_nthreads::Int = Base.Threads.nthreads(), 
    depletion_handling::Val{depletion_handling_enabled} = Val{false}(), 
    is_weighting_potential::Val{_is_weighting_potential} = Val{false}(),
    only2d::Val{only_2d} = Val{false}(),
)::Nothing where {T, depletion_handling_enabled, only_2d, _is_weighting_potential}
    outerloop!(pcs, use_nthreads, Val{true}(), depletion_handling, is_weighting_potential, only2d)
    apply_boundary_conditions!(pcs, Val{true}(), only2d)
    outerloop!(pcs, use_nthreads, Val{false}(), depletion_handling, is_weighting_potential, only2d)
    apply_boundary_conditions!(pcs, Val{false}(), only2d)
    nothing
end