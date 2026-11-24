using Plots
using Unitful
using SolidStateDetectors
T = Float64

mm = T(1/1000)
det_rin = 1*mm
det_z = det_r = 10*mm
z_draw = det_z/2
pn_r = 8.957282*mm # this one was calculated by searching the zero impurity point (displayed in the following section)

sim = Simulation{T}(SSD_examples[:TrueCoaxial])
cfn = SSD_examples[:TrueCoaxial]
print(open(f -> read(f, String), cfn))
plot(sim.detector, xunit = u"mm", yunit = u"mm", zunit = u"mm")
savefig("tutorial_det_dl.pdf") # hide

r_list = (0*mm):(0.01*mm):det_r
imp_list = T[]
for r in r_list
    pt = CylindricalPoint{T}(r, 0, z_draw)
    push!(imp_list, SolidStateDetectors.get_impurity_density(sim.detector.semiconductor.impurity_density_model, pt))
end
imp_list = imp_list ./ 1e6 # in cm^-3
plot(r_list/mm, imp_list, xlabel = "r [mm]", ylabel = "Impurity density [cm\$^{-3}\$]", label = "",
    color = :darkblue, lw = 2, grid = :on, xlims = (0, 10), ylims = (-2e10, 1e11))
vline!([pn_r/mm], lw = 2, ls = :dash, color = :darkred, label = "PN junction boundary")
savefig("tutorial_imp_dl.pdf") # hide

using SolidStateDetectors: Electron, Hole
cdm = sim.detector.semiconductor.charge_drift_model
depth_list = 0:(0.01*mm):(det_r-pn_r)
hole_mobility_list = []
electron_mobility_list = []
for depth in depth_list
    r = det_r-depth
    pt = CartesianPoint{T}(r, 0, z_draw)
    push!(hole_mobility_list, SolidStateDetectors.calculate_mobility(cdm, pt, Hole))
    push!(electron_mobility_list, SolidStateDetectors.calculate_mobility(cdm, pt, Electron))
end
plot(depth_list/mm, hole_mobility_list, label = "Hole", lw = 4)
plot!(depth_list/mm, electron_mobility_list, label = "Electron", lw = 4)
plot!(ylabel = "Mobility [cm\$^2\$/V/s]", xlabel = "Depth to surface [mm]",
    legend = :topleft, frame = :box, grid = :on, minorgrid = :on, xticks = 0:0.2:1.2, yticks = 0:0.5:4.5, ylims = [0, 4.5])
savefig("tutorial_mob_dl.pdf") # hide

calculate_electric_potential!(sim, max_n_iterations = 10, grid = Grid(sim), verbose = false, depletion_handling = true)
g = sim.electric_potential.grid
ax1, ax2, ax3 = g.axes
bulk_tick_dis = 0.05*mm
dl_tick_dis = 0.01*mm
user_additional_ticks_ax1 = sort(vcat(ax1.interval.left:bulk_tick_dis:pn_r, pn_r:dl_tick_dis:ax1.interval.right))
user_ax1 = typeof(ax1)(ax1.interval, user_additional_ticks_ax1)
user_g = typeof(g)((user_ax1, ax2, ax3))
calculate_electric_potential!(sim, refinement_limits = 0.1, use_nthreads = 8, grid = user_g, depletion_handling = true)
calculate_electric_field!(sim)
calculate_weighting_potential!(sim, 1, use_nthreads = 8, depletion_handling = true)
calculate_weighting_potential!(sim, 2, use_nthreads = 8, depletion_handling = true);
plot(
    begin
        imp = plot(sim.imp_scale, φ = 0, xunit = u"mm", yunit = u"mm", title = "impurity scale")
        vline!([pn_r/mm], lw = 2, ls = :dash, color = :darkred, label = "PN junction boundary", legendfontsize = 6)
    end,
    begin
        plot(sim.point_types, φ = 0, xunit = u"mm", yunit = u"mm", title = "point types")
        vline!([pn_r/mm], lw = 2, ls = :dash, color = :darkred, label = "PN junction boundary", legendfontsize = 6)
    end,
    begin
        plot(sim.electric_potential, title = "electric potential")
        vline!([pn_r/mm], lw = 2, ls = :dash, color = :darkred, label = "PN junction boundary", legendfontsize = 6)
    end,
    begin
        plot(sim.electric_field, title = "electric field", clims = (0, 100*2000))
        vline!([pn_r/mm], lw = 2, ls = :dash, color = :darkred, label = "PN junction boundary", legendfontsize = 6)
    end,
    size = (800, 600), layout = (2, 2),
)
savefig("tutorial_ef_dl.pdf") # hide

depth_list = (0.1*mm):(0.1*mm):(det_r-pn_r)
pulse_list = []
totTime = 5 # us
totEnergy = 1 # 1 keV --> simulating ~339 carrier pairs
max_nsteps = Int(totTime * 1000)
N = Int(totEnergy*1000÷2.95)
for depth in depth_list
    r = det_r-depth
    energy_depos = fill(T(2.95/1000), N) * u"keV"
    starting_positions = repeat([CartesianPoint{T}(r, deg2rad(0), z_draw)], N)
    evt = Event(starting_positions, energy_depos);
    simulate!(evt, sim, Δt = 1u"ns", max_nsteps = max_nsteps, diffusion = true, end_drift_when_no_field = false, self_repulsion = false)
    charge = ustrip(evt.waveforms[1].signal)
    push!(pulse_list, charge)
end
pulse_plot = plot()
eff_list = []
for (i, depth) in enumerate(depth_list)
    depth = round(depth/mm, digits = 1)
    plot!(pulse_list[i], label = "Depth: $(depth) mm", lw = 2, xlabel = "Time [ns]", ylabel = "Amplitude [e]",
        legend = :topright, grid = :on, minorgrid = :on, frame = :box)
    push!(eff_list, maximum(pulse_list[i])/N)
end
cce_plot = plot(depth_list/mm, eff_list, xlabel = "Depth to surface [mm]", lw = 2, ylabel = "Charge collection efficiency",
    frame = :box, grid = :on, minorgrid = :on, xticks = 0:0.2:1.2, yticks = 0:0.1:1, dpi = 500, color = :black, label = "")
plot(pulse_plot, cce_plot, layout = (1, 2), size = (1000, 400), margin = 5Plots.mm)
savefig("tutorial_pulse_cce_dl.pdf") # hide
