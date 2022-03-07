#!/bin/bash
#=
exec julia-1.6 --color=yes --startup-file=no --threads=auto "${BASH_SOURCE[0]}" "$@"
=#
# Execute this script via: `bash test/benchmarks/potential_calculation.jl overwrite`

# On Julia 1.7.X there seems to be an issue with KernelAbstractions.jl 
# -> https://github.com/JuliaGPU/KernelAbstractions.jl/issues/290
# such that the calculation hangs up sometimes.

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
fn_SOR_update_benchmark_df = joinpath(ENV["HOME"], "SSD_SOR_update_benchmark.jld2")

N_MAX_THREADS = Base.Threads.nthreads()

T = Float32

sim_cyl = Simulation{T}(SSD_examples[:BEGe])
sim_car = Simulation{T}(SSD_examples[:CGD])

hw_backends = (
    # (label = "CPU_1",  use_nthreads =   1, DAT = Array),
    # (label = "CPU_2",  use_nthreads =   2, DAT = Array),
    # (label = "CPU_4",  use_nthreads =   4, DAT = Array),
    (label = "CPU_8",  use_nthreads =   8, DAT = Array),
    # (label = "CPU_16", use_nthreads =  16, DAT = Array),
    # (label = "CPU_32", use_nthreads =  32, DAT = Array),
    # (label = "CPU_64", use_nthreads =  64, DAT = Array),
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
            verbose = true,
        )
        calculate_electric_potential!(sim; potential_calc_settings... )
        only2d = Val(size(sim.electric_potential, 2) == 1)

        for backend in hw_backends
            device_array_type = backend.DAT
            pcs = adapt(device_array_type, PotentialCalculationSetup(
                sim.detector, sim.electric_potential.grid, sim.medium, sim.electric_potential.data, sim.imp_scale.data, 
            ));
            ndrange = size(pcs.potential)[1:3] .- 2
            dev = get_device(pcs.potential)
            via_KernelAbstractions = false
            kernel = get_sor_kernel(S, dev, Val(via_KernelAbstractions))
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

display(df_pot_calc)
display(df_SOR_update)


# Plot the benchmarks

using Plots
pyplot()

function plot_benchmark(df::DataFrame; kwargs...)
    ys = df.time
    pos = eachindex(ys)
    plot(pos, ys,  
        xticks = (pos, df.backend),
        ylabel = "time [s]",
        st = :bar, 
        # ylims = (0.01, 1.1),
        # fillrange=1e-10,
        yscale = :identity,
        legend = false, 
        title = "$(df.S[1]) - Grid size: $(df.grid_size[1])";
        kwargs...
    )
end

plot(
    plot_benchmark(filter(row -> row.S == Cartesian,   df_pot_calc)),
    plot_benchmark(filter(row -> row.S == Cylindrical, df_pot_calc)),
    plot_benchmark(filter(row -> row.S == Cartesian,   df_SOR_update)),
    plot_benchmark(filter(row -> row.S == Cylindrical, df_SOR_update)),
    size = (1000, 1000),
    layout = (2, 2),
)

