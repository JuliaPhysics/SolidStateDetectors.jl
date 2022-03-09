using Plots; pyplot()
using DataFrames, JLD2
using SolidStateDetectors
using SolidStateDetectors: Cartesian, Cylindrical

fn_SOR_update_benchmark_df = joinpath(ENV["HOME"], "SSD_SOR_update_benchmark.jld2")
fn_pot_calc_benchmark_df = joinpath(ENV["HOME"], "SSD_potential_calculation_benchmark.jld2")

df_SOR_update = load(fn_SOR_update_benchmark_df, "df_SOR_update")
df_pot_calc = load(fn_pot_calc_benchmark_df, "df_SOR_update")

function plot_benchmark(df::DataFrame; kwargs...)
    ys = df.time 
    ys ./= minimum(ys)
    pos = eachindex(ys) 
    plot(pos, ys,  
        xticks = (pos, df.backend),
        ylabel = "Time [a.u.]",
        st = :bar, 
        # ylims = (0.001, max(ys)),
        # fillrange=1e-10,
        yscale = :identity,
        # yscale = :log10,
        legend = false, 
        title = "$(df.S[1]) - Grid size: $(df.grid_size[1])";
        kwargs...
    )
end

plot(
    plot_benchmark(filter(row -> row.S == Cartesian,   df_SOR_update)),
    plot_benchmark(filter(row -> row.S == Cylindrical, df_SOR_update)),
    size = (1200, 500),
    layout = (1, 2),
    guidefontsize = 16,
    tickfontsize = 12,
    titlefontsize = 18
)
savefig("/home/iwsatlas1/lhauert/Documents/sor_update_benchmark.pdf")
savefig("/home/iwsatlas1/lhauert/Documents/sor_update_benchmark.png")


plot(
    plot_benchmark(filter(row -> row.S == Cartesian,   df_SOR_update)),
    plot_benchmark(filter(row -> row.S == Cylindrical, df_SOR_update)),
    plot_benchmark(filter(row -> row.S == Cartesian,   df_pot_calc)),
    plot_benchmark(filter(row -> row.S == Cylindrical, df_pot_calc)),
    size = (1000, 1000),
    layout = (2, 2),
)


