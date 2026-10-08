# # Example 1: Inverted Coax Detector 

using Plots; pyplot(fmt = :png);
using SolidStateDetectors
using Unitful

T = Float32
simulation = Simulation{T}(SSD_examples[:InvertedCoax])

plot(simulation.detector, size = (700, 700))

# One can also have a look at how the initial conditions look like on the grid (its starts with a very coarse grid):

apply_initial_state!(simulation, ElectricPotential) # optional 
plot(
    plot(simulation.electric_potential), # initial electric potential (boundary conditions)
    plot(simulation.point_types), # map of different point types: fixed point / inside or outside detector volume / depleted/undepleted
    plot(simulation.ρ), # charge density distribution
    plot(simulation.ϵ), # dielectric distribution
    layout = (1, 4), size = (1400, 700)
)

# Next, calculate the electric potential:
 
calculate_electric_potential!( simulation, 
                               max_refinements = 3)

plot(
    plot(simulation.electric_potential, φ = 20), # initial electric potential (boundary conditions)
    plot(simulation.point_types), # map of different point types: fixed point / inside or outside detector volume / depleted/undepleted
    plot(simulation.ρ), # charge density distribution
    plot(simulation.ϵ), # dielectric distribution
    layout = (1, 4), size = (1400, 700)
)

# SolidStateDetectors.jl supports active (i.e. depleted) volume calculation:


get_active_volume(simulation.point_types) # approximation (sum of the volume of cells marked as depleted)

 
 
# ## Partially depleted detectors
 
# SolidStateDetectors.jl can also calculate the electric potential of a partially depleted detector:
 
detector_undep = deepcopy(simulation.detector)
detector_undep.contacts[end].potential = 500; # V  <-- Bias Voltage of Mantle

simulation_undep = Simulation(detector_undep);

calculate_electric_potential!( simulation_undep, 
                               depletion_handling = true, 
                               convergence_limit=1e-6,
                               max_refinements = 3, 
                               verbose = false)

plot(
    plot(simulation_undep.electric_potential), 
    plot(simulation_undep.point_types), 
    layout = (1, 2), size = (800, 700)
)

# Compare both volumes:

println("Depleted:   ", get_active_volume(simulation.point_types))
println("Undepleted: ", get_active_volume(simulation_undep.point_types));


# ## Electric field calculation

# Calculate the electric field of the fully depleted detector, given the already calculated electric potential:


calculate_electric_field!(simulation, n_points_in_φ = 72)

# plot_electric_field(simulation, size = (350, 500))

# ## Drift field calculation

# Given the electric field and a charge drift model, calculate drift fields for electrons and holes. Precalculating the drift fields saves time during charge drift simulation:

# Any drift field model can be used for the calculation of the electric field. If no model is explicitely given, the Bruyneel model from the Agata Data Library (ADL) is used. Other configurations are saved in their JSON configuration files and can be found under:

# `<package_directory>/src/ChargeDriftModels/ADL/<config_filename>.json.`

# Set the charge drift model of the simulation:

charge_drift_model = ADLChargeDriftModel()
set_charge_drift_model!(simulation, charge_drift_model) 

 
# And apply the charge drift model to the electric field:
 
calculate_drift_fields!(simulation)

# Now, let's create an "random" (multiside) event:

starting_positions = [ CylindricalPoint{T}( 0.020, deg2rad(10), 0.015 ), 
                       CylindricalPoint{T}( 0.015, deg2rad(20), 0.045 ), 
                       CylindricalPoint{T}( 0.025, deg2rad(30), 0.025 ) ]
energy_depos = T[1460, 609, 1000] * u"keV" # are needed later in the signal generation

event = Event(starting_positions, energy_depos);

time_step = 5u"ns"
drift_charges!(event, simulation, Δt = time_step)

plot(simulation.detector, size = (700, 700))
plot!(event.drift_paths)

# ## Weighting potential calculation

# We need weighting potentials to simulate the detector charge signal induced by drifting charges. We'll calculate the weighting potential for the point contact and the outer shell of the detector:

for contact in simulation.detector.contacts
    calculate_weighting_potential!(simulation, contact.id, max_refinements = 3, n_points_in_φ = 2, verbose = false)
end

plot(  
    plot(simulation.weighting_potentials[1]),
    plot(simulation.weighting_potentials[2]),
    size = (900, 700)
)

# ## Detector waveform generation

# ### Single-event simulation
 
# Given an interaction at an arbitrary point in the detector, we can now simulate the charge drift and the resulting detector charge signals (e.g. at the point contact):

simulate!(event, simulation) # drift_charges + signal generation of all channels

p_pc_signal = plot( event.waveforms[1], lw = 1.5, xlims = (0, 2000), xlabel = "Time / ns", 
                    legend = false, tickfontsize = 12, ylabel = "Energy / eV", guidefontsize = 14)
