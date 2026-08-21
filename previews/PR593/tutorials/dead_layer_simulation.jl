using Plots
using Unitful
using SolidStateDetectors
T = Float64;

det_rin = 1u"mm"
det_z = det_r = 10u"mm"
z_draw = det_z/2
pn_r = 8.957282u"mm" # this one was calculated by searching the zero impurity point (displayed in the following section)

sim = Simulation{T}(SSD_examples[:TrueCoaxial])
cfn = SSD_examples[:TrueCoaxial]
print(open(f -> read(f, String), cfn))
plot(sim.detector, xunit = u"mm", yunit = u"mm", zunit = u"mm")
savefig("tutorial_det_dl.pdf") # hide

r_list = 0u"mm":0.01u"mm":det_r
imp_list = map(r -> let pt::CylindricalPoint{T} = CylindricalPoint(r, 0u"°", z_draw)
    SolidStateDetectors.get_impurity_density(sim.detector.semiconductor.impurity_density_model, pt) * 1e-6u"cm^-3"
end, r_list)
plot(r_list, imp_list, xlabel = "r / mm", ylabel = "Impurity density / cm\$^{-3}\$", unitformat = :nounit, label = "",
    color = :darkblue, lw = 2, grid = :on, xlims = (0, 10), ylims = (-2e10, 1e11))
vline!([pn_r], lw = 2, ls = :dash, color = :darkred, label = "PN junction boundary")
savefig("tutorial_imp_dl.pdf") # hide

using SolidStateDetectors: Electron, Hole
cdm = sim.detector.semiconductor.charge_drift_model
depth_list = 0u"mm":0.01u"mm":(det_r-pn_r)
mobility_list = map(depth -> let pt::CartesianPoint{T} = CartesianPoint(det_r - depth, 0, z_draw)
    (µe = SolidStateDetectors.calculate_mobility(cdm, pt, Hole) * 10000u"cm^2/(V*s)",
     µh = SolidStateDetectors.calculate_mobility(cdm, pt, Electron) * 10000u"cm^2/(V*s)")
end, depth_list)
plot(depth_list,  getfield.(mobility_list, :µh), label = "Hole", lw = 4)
plot!(depth_list, getfield.(mobility_list, :µe), label = "Electron", lw = 4)
plot!(xlabel = "Depth to surface / mm", ylabel = "Mobility / cm\$^2\$/Vs", unitformat = :nounit, legend = :topleft, xlims = (0u"mm", det_r-pn_r))
savefig("tutorial_mob_dl.pdf") # hide

calculate_electric_potential!(sim, max_n_iterations = 10, grid = Grid(sim), verbose = false, depletion_handling = true)
g = sim.electric_potential.grid
ax1, ax2, ax3 = g.axes
bulk_tick_dis = 0.05u"mm"
dl_tick_dis   = 0.01u"mm"
user_additional_ticks_ax1 = sort(vcat(ax1.interval.left*u"m":bulk_tick_dis:pn_r, pn_r:dl_tick_dis:ax1.interval.right*u"m"))
user_ax1 = typeof(ax1)(ax1.interval, SolidStateDetectors.to_internal_units.(user_additional_ticks_ax1))
user_g = typeof(g)((user_ax1, ax2, ax3))
calculate_electric_potential!(sim, refinement_limits = 0.1, grid = user_g, depletion_handling = true)
calculate_electric_field!(sim)
calculate_weighting_potential!(sim, 1, depletion_handling = true)
calculate_weighting_potential!(sim, 2, depletion_handling = true);
plot(
    begin
        imp = plot(sim.imp_scale, φ = 0, xunit = u"mm", yunit = u"mm", title = "impurity scale")
        vline!([pn_r], lw = 2, ls = :dash, color = :darkred, label = "PN junction boundary", legendfontsize = 6)
    end,
    begin
        plot(sim.point_types, φ = 0, xunit = u"mm", yunit = u"mm", title = "point types")
        vline!([pn_r], lw = 2, ls = :dash, color = :darkred, label = "PN junction boundary", legendfontsize = 6)
    end,
    begin
        plot(sim.electric_potential, xunit = u"mm", yunit = u"mm", title = "electric potential")
        vline!([pn_r], lw = 2, ls = :dash, color = :darkred, label = "PN junction boundary", legendfontsize = 6)
    end,
    begin
        plot(sim.electric_field, xunit = u"mm", yunit = u"mm", title = "electric field", clims = (0, 100*2000))
        vline!([pn_r], lw = 2, ls = :dash, color = :darkred, label = "PN junction boundary", legendfontsize = 6)
    end,
    size = (800, 600), layout = (2, 2),
)
savefig("tutorial_ef_dl.pdf") # hide

depth_list = 0.1u"mm":0.1u"mm":(det_r-pn_r)
totTime = 5u"µs"
totEnergy = 1u"keV" # --> simulating ~339 carrier pairs
N = Int(totEnergy÷2.95u"eV")

pulse_plot = plot()
eff_list = map(depth -> begin
    r = det_r-depth
    energy_depos = fill(2.95u"eV", N)
    starting_positions = repeat([CartesianPoint(r, 0, z_draw)], N)
    evt = Event(starting_positions, energy_depos)
    simulate!(evt, sim, Δt = 1u"ns", max_nsteps = round(Int, totTime / 1u"ns"), diffusion = true, end_drift_when_no_field = false, self_repulsion = false)
    pulse = evt.waveforms[1]
    plot!(pulse_plot, pulse, label = "Depth: $(round(typeof(depth), depth, digits = 1))", lw = 2, yunit = u"eV/V")
    maximum(pulse.signal)/N
end, depth_list)
plot!(pulse_plot, legend = :topright, xlabel = "Time / ns", ylabel = "Amplitude / e", unitformat = :nounit)
cce_plot = plot(depth_list, eff_list, xlabel = "Depth to surface / mm", ylabel = "Charge collection efficiency", lw = 2, color = :black, label = "", unitformat = :nounit)
plot(pulse_plot, cce_plot, layout = (1, 2), size = (1000, 400), margin = 5Plots.mm)
savefig("tutorial_pulse_cce_dl.pdf") # hide
