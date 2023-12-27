#!/bin/bash
#=
exec julia --color=yes --startup-file=no --threads=auto "${BASH_SOURCE[0]}" "$@"
=#
# Execute this script via: `bash test/benchmarks/benchmark_SOR_update_method.jl overwrite`

recalculate = in("overwrite", ARGS)

using CUDAKernels, CUDA
using SolidStateDetectors
using DataFrames, JLD2
using BenchmarkTools
using KernelAbstractions

using SolidStateDetectors: 
    ConstructiveSolidGeometry.AbstractCoordinateSystem,
    Cartesian, Cylindrical, PotentialCalculationSetup,
    _guess_optimal_number_of_threads_for_SOR,
    get_sor_kernel, update!, _ka_get_backend

import Adapt

fn_pot_calc_benchmark_df = joinpath(ENV["HOME"], "SSD_potential_calculation_benchmark.jld2")
fn_SOR_update_benchmark_df = joinpath(ENV["HOME"], "SSD_SOR_update_benchmark.jld2")

N_MAX_THREADS = Base.Threads.nthreads()

T = Float32

sim_cyl = Simulation{T}(SSD_examples[:BEGe])
sim_car = Simulation{T}(SSD_examples[:CGD])

hw_backends = (
    (label = "CPU_1",  use_nthreads =   1, DAT = Array),
    (label = "CPU_2",  use_nthreads =   2, DAT = Array),
    (label = "CPU_4",  use_nthreads =   4, DAT = Array),
    (label = "CPU_8",  use_nthreads =   8, DAT = Array),
    (label = "CPU_16", use_nthreads =  16, DAT = Array),
    (label = "CPU_32", use_nthreads =  32, DAT = Array),
    (label = "GPU",    use_nthreads =   1, DAT = CuArray),
)

refinement_limits = [0.2, 0.1, 0.05, 0.025, 0.01]
depletion_handling = Val(true)
is_weighting_potential = Val(false)

df_SOR_update = if !isfile(fn_SOR_update_benchmark_df) || recalculate
    df = DataFrame()
    for sim in (sim_cyl, sim_car)
        S = SolidStateDetectors.get_coordinate_system(sim)
        potential_calc_settings = (
            convergence_limit = 1e-5,
            refinement_limits = refinement_limits,
            depletion_handling = true,
            use_nthreads = 16,
            verbose = true,
        )
        calculate_electric_potential!(sim; potential_calc_settings... )
        only2d = Val(size(sim.electric_potential, 2) == 1)

        for backend in hw_backends
            device_array_type = backend.DAT
            pcs = Adapt.adapt(device_array_type, PotentialCalculationSetup(
                sim.detector, sim.electric_potential.grid, sim.medium, sim.electric_potential.data, sim.imp_scale.data, 
            ));
            ndrange = size(pcs.potential)[1:3] .- 2
            backend = _ka_get_backend(pcs.potential)
            via_KernelAbstractions = false
            kernel = get_sor_kernel(S, backend, Val(via_KernelAbstractions))
            use_nthreads = backend.use_nthreads
            update!(pcs, kernel, ndrange; use_nthreads, depletion_handling, is_weighting_potential, only2d)
            bm = @benchmark update!($pcs, $kernel, $ndrange; 
                use_nthreads = $use_nthreads, 
                depletion_handling = $depletion_handling, 
                is_weighting_potential = $is_weighting_potential, 
                only2d = $only2d
            )
            time_μ = BenchmarkTools.mean(bm.times) * 1e-9 # now in seconds
            time_σ = BenchmarkTools.std(bm.times)  * 1e-9 # now in seconds
            push!(df, 
                (
                    backend = backend.label,
                    S = S,
                    time = time_μ,
                    time_σ = time_σ,
                    grid_size = size(sim.electric_potential),
                    potential_calc_settings...
                )
            )
            display(df)
        end
    end
    save(fn_SOR_update_benchmark_df, Dict("df_SOR_update" => df))
    df
else
    load(fn_SOR_update_benchmark_df, "df_SOR_update")
end