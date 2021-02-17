using Plots; pyplot(fmt = :png);
using SolidStateDetectors
using Unitful

T = Float32
simulation = Simulation{T}(SSD_examples[:InvertedCoax])

plot(simulation.detector, size = (700, 700))

apply_initial_state!(simulation, ElectricPotential) # optional
plot(
    plot(simulation.electric_potential), # initial electric potential (boundary conditions)
    plot(simulation.point_types), # map of different point types: fixed point / inside or outside detector volume / depleted/undepleted
    plot(simulation.ρ), # charge density distribution
    plot(simulation.ϵ), # dielectric distribution
    layout = (1, 4), size = (1400, 700)
)

calculate_electric_potential!( simulation,
                               max_refinements = 4)

plot(
    plot(simulation.electric_potential, φ = 20), # initial electric potential (boundary conditions)
    plot(simulation.point_types), # map of different point types: fixed point / inside or outside detector volume / depleted/undepleted
    plot(simulation.ρ), # charge density distribution
    plot(simulation.ϵ), # dielectric distribution
    layout = (1, 4), size = (1400, 700)
)

get_active_volume(simulation.point_types) # approximation (sum of the volume of cells marked as depleted)

detector_undep = deepcopy(simulation.detector)
detector_undep.contacts[end].potential = 500; # V  <-- Bias Voltage of Mantle

simulation_undep = Simulation(detector_undep);

calculate_electric_potential!( simulation_undep,
                               depletion_handling = true,
                               max_refinements = 4,
                               verbose = false)

plot(
    plot(simulation_undep.electric_potential),
    plot(simulation_undep.point_types),
    layout = (1, 2), size = (800, 700)
)

println("Depleted:   ", get_active_volume(simulation.point_types))
println("Undepleted: ", get_active_volume(simulation_undep.point_types));

calculate_electric_field!(simulation, n_points_in_φ = 72)

plot_electric_field(simulation, size = (350, 500))

apply_charge_drift_model!(simulation)

starting_positions = [ CylindricalPoint{T}( 0.020, deg2rad(10), 0.015 ),
                       CylindricalPoint{T}( 0.015, deg2rad(20), 0.045 ),
                       CylindricalPoint{T}( 0.025, deg2rad(30), 0.025 ) ]
energy_depos = T[1460, 609, 1000] * u"keV" # are needed later in the signal generation

event = Event(starting_positions, energy_depos);

time_step = 5u"ns"
drift_charges!(event, simulation, Δt = time_step)

plot(simulation.detector, size = (700, 700))
plot!(event.drift_paths)

for contact in simulation.detector.contacts
    calculate_weighting_potential!(simulation, contact.id, max_refinements = 4, n_points_in_φ = 2, verbose = false)
end

plot(
    plot(simulation.weighting_potentials[1]),
    plot(simulation.weighting_potentials[2]),
    size = (900, 700)
)

simulate!(event, simulation) # drift_charges + signal generation of all channels

p_pc_signal = plot( event.waveforms[1], lw = 1.5, xlims = (0, 2000), xlabel = "Time / ns",
                    legend = false, tickfontsize = 12, ylabel = "Energy / eV", guidefontsize = 14)

