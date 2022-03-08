#!/bin/bash
#=
exec julia --color=yes --startup-file=no --threads=auto "${BASH_SOURCE[0]}" "$@"
=#
# Execute this script via: `bash test/benchmarks/benchmark_potential_calculation.jl overwrite`

recalculate = in("overwrite", ARGS)

using CUDAKernels, CUDA
using SolidStateDetectors
using DataFrames, JLD2
using BenchmarkTools
using KernelAbstractions
using KernelAbstractions: get_device

using SolidStateDetectors: 
    ConstructiveSolidGeometry.AbstractCoordinateSystem,
    Cartesian, Cylindrical, PotentialCalculationSetup,
    _guess_optimal_number_of_threads_for_SOR, adapt,
    get_sor_kernel, update!

fn_pot_calc_benchmark_df = joinpath(ENV["HOME"], "SSD_potential_calculation_benchmark.jld2")

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
max_n_iterations = 10_000

df_pot_calc = if !isfile(fn_pot_calc_benchmark_df) || recalculate
    df = DataFrame()
    for sim in (sim_cyl, sim_car)
        S = SolidStateDetectors.get_coordinate_system(sim)
        for backend in hw_backends
            potential_calc_settings = (
                use_nthreads = backend.use_nthreads,
                convergence_limit = 1e-6,
                refinement_limits = refinement_limits,
                device_array_type = backend.DAT,
                max_n_iterations = max_n_iterations,
                depletion_handling = true,
                verbose = true,
            )
            calculate_electric_potential!(sim; # compile time
                max_n_iterations = 40, 
                n_iterations_between_checks = 10,
                potential_calc_settings...
            )                
            t = (@elapsed calculate_electric_potential!(sim; potential_calc_settings... ))
            push!(df, 
                (
                    backend = backend.label,
                    S = S,
                    time = t,
                    grid_size = size(sim.electric_potential),
                    potential_calc_settings...
                )
            )
            display(df)
        end
    end
    save(fn_pot_calc_benchmark_df, Dict("df_pot_calc" => df))
    df
else
    load(fn_pot_calc_benchmark_df, "df_pot_calc")
end
