using SolidStateDetectors
using LegendHDF5IO
using Unitful
using JLD2
using Plots
using BenchmarkTools

sim_file = "$(@__DIR__)/../../ssd-sims.h5"
sim = if ispath(sim_file)
    ssd_read(sim_file, Simulation)
else
    @info "loading simulation didn't work"
end 

T = Float32

starting_positions = [CartesianPoint{T}(-0.02, 0.015, 0.04), 
    CartesianPoint{T}(0.015, -0.012, 0.02), 
    CartesianPoint{T}(0.01, -0.025, 0.01)]
energy_depos = T[1460, 609, 1000] * u"keV" # are needed later in the signal generation


evt = Event(starting_positions, energy_depos, 4)

time_step = 5u"ns"

drift_charges!(evt, sim, Δt = time_step; self_repulsion = true)

new_drift_paths = evt.drift_paths

# run original simulation chain to compare results
# original_drift_paths = JLD2.load("$(@__DIR__)/../../ssd_simulation_chain_original.jld2", "drift_paths")


plot(sim.detector, xunit = u"mm", yunit = u"mm", zunit = u"mm", size = (700, 700))
plot!(original_drift_paths, color = :red, linewidth = 2)
plot!(new_drift_paths, color = :blue, linewidth = 2)

for i in eachindex(new_drift_paths)
    a = findall(new_drift_paths[i].e_path != original_drift_paths[i].e_path)
    if a != []
        @show i
        @show a
    else nothing
    end
    b = findall(new_drift_paths[i].h_path != original_drift_paths[i].h_path)
    if b != []
        @show i
        @show b
    else nothing
    end
    c = findall(new_drift_paths[i].timestamps_e != original_drift_paths[i].timestamps_e)
    if c != []
        @show i
        @show c
    else nothing
    end
    d = findall(new_drift_paths[i].timestamps_h != original_drift_paths[i].timestamps_h)
    if d != []
        @show i
        @show  d
    else nothing
    end
end


# new drift paths does not have the same length as original drift paths
tmp1 = new_drift_paths[1].e_path
tmp2 = original_drift_paths[1].e_path

length(new_drift_paths[1].e_path)
length(original_drift_paths[1].e_path) # should be true