# # Example 1: Inverted Coax Detector

using Plots
using SolidStateDetectors
using Unitful

T = Float32
simulation = Simulation{T}(SSD_examples[:InvertedCoax])

plot(simulation.detector)
#jl savefig("tutorial_det.pdf") # hide
#md savefig("tutorial_det.pdf") # hide
#md savefig("tutorial_det.svg"); nothing # hide
#md # [![tutorial_det](tutorial_det.svg)](tutorial_det.pdf)

# One can also have a look at how the initial conditions look like on the grid (its starts with a very coarse grid):

apply_initial_state!(simulation, ElectricPotential) # optional
plot(
    plot(simulation.electric_potential), # initial electric potential (boundary conditions)
    plot(simulation.point_types), # map of different point types: fixed point / inside or outside detector volume / depleted/undepleted
    plot(simulation.q_eff_imp), # charge density distribution
    plot(simulation.ϵ_r), # dielectric distribution
    layout = (1, 4), size = (1600, 500)
)
#jl savefig("tutorial_initial_condition.pdf") # hide
#md savefig("tutorial_initial_condition.pdf") # hide
#md savefig("tutorial_initial_condition.svg"); nothing # hide
#md # [![tutorial_initial_condition](tutorial_initial_condition.svg)](tutorial_initial_condition.pdf)


# Next, calculate the electric potential:

calculate_electric_potential!( simulation,
                               refinement_limits = [0.2, 0.1, 0.05, 0.01])

plot(
    plot(simulation.electric_potential, φ = 20), # initial electric potential (boundary conditions)
    plot(simulation.point_types), # map of different point types: fixed point / inside or outside detector volume / depleted/undepleted
    plot(simulation.q_eff_imp), # charge density distribution
    plot(simulation.ϵ_r), # dielectric distribution
    layout = (1, 4), size = (1600, 500)
)
#jl savefig("tutorial_calculated_potential.pdf") # hide
#md savefig("tutorial_calculated_potential.pdf") # hide
#md savefig("tutorial_calculated_potential.svg"); nothing # hide
#md # [![tutorial_calculated_potential](tutorial_calculated_potential.svg)](tutorial_calculated_potential.pdf)

# SolidStateDetectors.jl supports active (i.e. depleted) volume calculation:


get_active_volume(simulation.point_types) # approximation (sum of the volume of cells marked as depleted)



# ## Partially depleted detectors

# SolidStateDetectors.jl can also calculate the electric potential of a partially depleted detector:

simulation_undep = deepcopy(simulation)
simulation_undep.detector = SolidStateDetector(simulation_undep.detector, contact_id = 2, contact_potential = 500); # V  <-- Bias Voltage of Mantle

calculate_electric_potential!( simulation_undep,
                               depletion_handling = true,
                               convergence_limit=1e-6,
                               refinement_limits = [0.2, 0.1, 0.05, 0.01],
                               verbose = false)


plot(
    plot(simulation_undep.electric_potential),
    plot(simulation_undep.point_types),
    layout = (1, 2), size = (800, 700)
)
#jl savefig("tutorial_calculated_potential_undep.pdf") # hide
#md savefig("tutorial_calculated_potential_undep.pdf") # hide
#md savefig("tutorial_calculated_potential_undep.svg"); nothing # hide
#md # [![tutorial_calculated_potential_undep](tutorial_calculated_potential_undep.svg)](tutorial_calculated_potential_undep.pdf)

# Compare both volumes:

println("Depleted:   ", get_active_volume(simulation.point_types))
println("Undepleted: ", get_active_volume(simulation_undep.point_types));


# ## Electric field calculation

# Calculate the electric field of the fully depleted detector, given the already calculated electric potential:


calculate_electric_field!(simulation, n_points_in_φ = 72)

plot(simulation.electric_field, full_det = true, φ = 0.0, size = (500, 500))
plot_electric_fieldlines!(simulation, full_det = true, φ = 0.0)
#jl savefig("tutorial_electric_field.pdf") # hide
#md savefig("tutorial_electric_field.pdf") # hide
#md savefig("tutorial_electric_field.svg"); nothing # hide
#md # [![tutorial_electric_field](tutorial_electric_field.svg)](tutorial_electric_field.pdf)


# ## Drift field calculation

# Given the electric field and a charge drift model, calculate drift fields for electrons and holes. Precalculating the drift fields saves time during charge drift simulation:

# Any drift field model can be used for the calculation of the electric field. If no model is explicitely given, the Bruyneel model from the Agata Data Library (ADL) is used. Other configurations are saved in their configuration files and can be found under:

# `<package_directory>/examples/example_config_files/ADLChargeDriftModel/<config_filename>.yaml.`

# Set the charge drift model of the simulation:

charge_drift_model = ADLChargeDriftModel()
simulation.detector = SolidStateDetector(simulation.detector, charge_drift_model)


# And apply the charge drift model to the electric field:

calculate_drift_fields!(simulation)

# Now, let's create an "random" (multiside) event:

starting_positions = [ CylindricalPoint{T}( 0.020, deg2rad(10), 0.015 ),
                       CylindricalPoint{T}( 0.015, deg2rad(20), 0.045 ),
                       CylindricalPoint{T}( 0.022, deg2rad(35), 0.025 ) ]
energy_depos = T[1460, 609, 1000] * u"keV" # are needed later in the signal generation

event = Event(starting_positions, energy_depos);

time_step = 5u"ns"
drift_charges!(event, simulation, Δt = time_step)

plot(simulation.detector, size = (700, 700))
plot!(event.drift_paths)
#jl savefig("tutorial_drift_paths.pdf") # hide
#md savefig("tutorial_drift_paths.pdf") # hide
#md savefig("tutorial_drift_paths.svg"); nothing # hide
#md # [![tutorial_drift_paths](tutorial_drift_paths.svg)](tutorial_drift_paths.pdf)


# ## Weighting potential calculation

# We need weighting potentials to simulate the detector charge signal induced by drifting charges. We'll calculate the weighting potential for the point contact and the outer shell of the detector:

for contact in simulation.detector.contacts
    calculate_weighting_potential!(simulation, contact.id, refinement_limits = [0.2, 0.1, 0.05, 0.01], n_points_in_φ = 2, verbose = false)
end

plot(
    plot(simulation.weighting_potentials[1]),
    plot(simulation.weighting_potentials[2]),
    size = (900, 700)
)
#jl savefig("tutorial_weighting_potentials.pdf") # hide
#md savefig("tutorial_weighting_potentials.pdf") # hide
#md savefig("tutorial_weighting_potentials.svg"); nothing # hide
#md # [![tutorial_weighting_potentials](tutorial_weighting_potentials.svg)](tutorial_weighting_potentials.pdf)


# ## Detector waveform generation

# ### Single-event simulation

# Given an interaction at an arbitrary point in the detector, we can now simulate the charge drift and the resulting detector charge signals (e.g. at the point contact):

simulate!(event, simulation) # drift_charges + signal generation of all channels

p_pc_signal = plot( event.waveforms[1], lw = 1.5, xlims = (0, 1100), xlabel = "Time / ns",
                    legend = false, tickfontsize = 12, ylabel = "Energy / eV", guidefontsize = 14)
#jl savefig("tutorial_waveforms.pdf") # hide
#md savefig("tutorial_waveforms.pdf") # hide
#md savefig("tutorial_waveforms.svg"); nothing # hide
#md # [![tutorial_waveforms](tutorial_waveforms.svg)](tutorial_waveforms.pdf)

