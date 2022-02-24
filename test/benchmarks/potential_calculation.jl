using CUDAKernels
using CUDAKernels.CUDA: CuArray
using SolidStateDetectors
using DataFrames, JLD2
using Unitful

using SolidStateDetectors: 
    ConstructiveSolidGeometry.AbstractCoordinateSystem,
    Cartesian, Cylindrical

fn_pot_calc_benchmark_df = joinpath(ENV["HOME"], "SSD_potential_calculation_benchmark.jld2")

N_MAX_THREADS = Base.Threads.nthreads()

T = Float32

sim_cyl = Simulation{T}(SSD_examples[:BEGe])
sim_car = Simulation{T}(SSD_examples[:CGD])

max_n_threads = [1, 2, 4, 8, 16, 32, 64, missing]
max_refinement_limits = [0.2, 0.1, 0.05]#, 0.025, 0.01, 0.005]

overwrite = false

df_pot_calc = if !isfile(fn_pot_calc_benchmark_df) || overwrite
    df = DataFrame()
    for sim in (sim_cyl, sim_car)
        S = SolidStateDetectors.get_coordinate_system(sim)
        for nt in max_n_threads
            device_array_type = ismissing(nt) ? CuArray : Array # GPU : CPU
            max_nt = ismissing(nt) ? 1000 : nt # not used for GPU as backend -> just set it to 1000
            refinement_limits = [0.2]
            for i_final_ref in eachindex(max_refinement_limits)  
                potential_calc_settings = (
                    use_nthreads = max_nt,
                    convergence_limit = 0,
                    refinement_limits = max_refinement_limits[1:i_final_ref],
                    depletion_handling = false,
                    max_n_iterations = 100,
                    verbose = true,
                )
                t = (@elapsed calculate_electric_potential!(sim; potential_calc_settings... ))
                push!(df, 
                    (
                        S = S,
                        time = t,
                        DAT = device_array_type,
                        potential_calc_settings...
                    )
                )
                display(df)
            end
        end
    end
    save(fn_pot_calc_benchmark_df, Dict("df_pot_calc" => df))
    df
else
    load(fn_pot_calc_benchmark_df, "df_pot_calc")
end

using Plots, StatsPlots

begin
    subdf = filter(row -> row.S == Cartesian, df_pot_calc)

    df_refs = groupby(subdf, :refinement_limits)
    
    n_cats = size(df_refs[1], 1)
    n_groups = length(df_refs)
    data = zeros(n_cats, n_groups)    
    for iG in 1:n_groups
        df = df_refs[iG]
        data[:, iG] = df[!, :time]
    end

    group_names = string.(repeat(
        [df.refinement_limits[1][end] for df in df_refs], 
        outer = size(data, 1)
    ))
    category_names = repeat(
        vcat(("CPU: " .* (string.(df_refs[1].use_nthreads[1:end-1]))), ["GPU"]), 
        inner = size(data, 2)
    )
    
    groupedbar(group_names, data, group = category_names, 
        title = "Potential Calculation Benchmark", bar_width = 0.8,
        xlabel = "Final ref. limit [1]", ylabel = "Time [s]",
        lw = 1, framestyle = :box
    )


    groupedbar(group_names, data, group = category_names, xlabel = "Final ref. limit [1]", ylabel = "Time [s]", 
            title = "Potential Calculation Benchmark", bar_width = 0.1,
            lw = 1, framestyle = :box)
end

