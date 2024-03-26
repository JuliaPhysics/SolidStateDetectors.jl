using Plots
using SolidStateDetectors
using Unitful

T = Float32
sim = Simulation{T}(SSD_examples[:InvertedCoax])

plot(sim.detector, size = (700, 700))
savefig("tutorial_det.pdf") # hide

apply_initial_state!(sim, ElectricPotential) # optional
plot(
    plot(sim.electric_potential), # initial electric potential (boundary conditions)
    plot(sim.point_types), # map of different point types: fixed point / inside or outside detector volume / depleted/undepleted
    plot(sim.q_eff_imp), # charge density distribution
    plot(sim.ϵ_r), # dielectric distribution
    layout = (1, 4), size = (1600, 500)
)
savefig("tutorial_initial_condition.pdf") # hide

calculate_electric_potential!( sim,
                               refinement_limits = [0.2, 0.1, 0.05, 0.01])

plot(
    plot(sim.electric_potential, φ = 20), # initial electric potential (boundary conditions)
    plot(sim.point_types), # map of different point types: fixed point / inside or outside detector volume / depleted/undepleted
    plot(sim.q_eff_imp), # charge density distribution
    plot(sim.ϵ_r), # dielectric distribution
    layout = (1, 4), size = (1600, 500)
)
savefig("tutorial_calculated_potential.pdf") # hide

is_depleted(sim.point_types)

get_active_volume(sim.point_types) # approximation (sum of the volume of cells marked as depleted)

sim_undep = deepcopy(sim)
sim_undep.detector = SolidStateDetector(sim_undep.detector, contact_id = 2, contact_potential = 500); # V  <-- Bias Voltage of Mantle

calculate_electric_potential!( sim_undep,
                               depletion_handling = true,
                               convergence_limit = 1e-6,
                               refinement_limits = [0.2, 0.1, 0.05, 0.01],
                               verbose = false)


plot(
    plot(sim_undep.electric_potential),
    plot(sim_undep.point_types),
    layout = (1, 2), size = (800, 700)
)
savefig("tutorial_calculated_potential_undep.pdf") # hide

is_depleted(sim_undep.point_types)

println("Depleted:   ", get_active_volume(sim.point_types))
println("Undepleted: ", get_active_volume(sim_undep.point_types));

calculate_electric_field!(sim, n_points_in_φ = 72)

plot(sim.electric_field, full_det = true, φ = 0.0, size = (700, 700))
plot_electric_fieldlines!(sim, full_det = true, φ = 0.0)
savefig("tutorial_electric_field.pdf") # hide

charge_drift_model = ADLChargeDriftModel()
sim.detector = SolidStateDetector(sim.detector, charge_drift_model)

starting_positions = [ CylindricalPoint{T}( 0.020, deg2rad(10), 0.015 ),
                       CylindricalPoint{T}( 0.015, deg2rad(20), 0.045 ),
                       CylindricalPoint{T}( 0.022, deg2rad(35), 0.025 ) ]
energy_depos = T[1460, 609, 1000] * u"keV" # are needed later in the signal generation

evt = Event(starting_positions, energy_depos);

time_step = 5u"ns"
drift_charges!(evt, sim, Δt = time_step)

plot(sim.detector, size = (700, 700))
plot!(evt.drift_paths)
savefig("tutorial_drift_paths.pdf") # hide

for contact in sim.detector.contacts
    calculate_weighting_potential!(sim, contact.id, refinement_limits = [0.2, 0.1, 0.05, 0.01], n_points_in_φ = 2, verbose = false)
end

plot(
    plot(sim.weighting_potentials[1]),
    plot(sim.weighting_potentials[2]),
    size = (900, 700)
)
savefig("tutorial_weighting_potentials.pdf") # hide

simulate!(evt, sim) # drift_charges + signal generation of all channels

p_pc_signal = plot( evt.waveforms[1], lw = 1.5, xlims = (0, 1100), xlabel = "Time / ns",
                    legend = false, tickfontsize = 12, ylabel = "Energy / eV", guidefontsize = 14)
savefig("tutorial_waveforms.pdf") # hide

contact_id = 1
plot_electron_and_hole_contribution(evt, sim, contact_id, xlims = (0, 1100), xlabel = "Time / ns",
                    legend = :topleft, tickfontsize = 12, ylabel = "Energy / eV", guidefontsize = 14)
savefig("tutorial_waveform_contributions.pdf") # hide

