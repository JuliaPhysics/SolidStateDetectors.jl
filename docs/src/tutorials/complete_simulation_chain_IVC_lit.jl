# # Simulation Chain: Inverted Coax Detector

using Plots
using SolidStateDetectors
using Unitful

T = Float32
sim = Simulation{T}(SSD_examples[:InvertedCoax])

plot(sim.detector, size = (700, 700))
#jl savefig("tutorial_det.pdf") # hide
#md savefig("tutorial_det.pdf") # hide
#md savefig("tutorial_det.svg"); nothing # hide
#md # [![tutorial_det](tutorial_det.svg)](tutorial_det.pdf)

# One can also have a look at how the initial conditions look like on the grid (its starts with a very coarse grid):

apply_initial_state!(sim, ElectricPotential) # optional
plot(
    plot(sim.electric_potential), # initial electric potential (boundary conditions)
    plot(sim.point_types), # map of different point types: fixed point / inside or outside detector volume / depleted/undepleted
    plot(sim.q_eff_imp), # charge density distribution
    plot(sim.ϵ_r), # dielectric distribution
    layout = (1, 4), size = (1600, 500)
)
#jl savefig("tutorial_initial_condition.pdf") # hide
#md savefig("tutorial_initial_condition.pdf") # hide
#md savefig("tutorial_initial_condition.svg"); nothing # hide
#md # [![tutorial_initial_condition](tutorial_initial_condition.svg)](tutorial_initial_condition.pdf)


# Next, calculate the electric potential:

calculate_electric_potential!( sim,
                               refinement_limits = [0.2, 0.1, 0.05, 0.01])

plot(
    plot(sim.electric_potential, φ = 20), # initial electric potential (boundary conditions)
    plot(sim.point_types), # map of different point types: fixed point / inside or outside detector volume / depleted/undepleted
    plot(sim.q_eff_imp), # charge density distribution
    plot(sim.ϵ_r), # dielectric distribution
    layout = (1, 4), size = (1600, 500)
)
#jl savefig("tutorial_calculated_potential.pdf") # hide
#md savefig("tutorial_calculated_potential.pdf") # hide
#md savefig("tutorial_calculated_potential.svg"); nothing # hide
#md # [![tutorial_calculated_potential](tutorial_calculated_potential.svg)](tutorial_calculated_potential.pdf)

# SolidStateDetectors.jl supports active (i.e. depleted) volume calculation:

is_depleted(sim.point_types)

# 

get_active_volume(sim.point_types) # approximation (sum of the volume of cells marked as depleted)



# ## Partially depleted detectors

# SolidStateDetectors.jl can also calculate the electric potential of a partially depleted detector:

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
    plot(sim_undep.imp_scale),
    layout = (1, 3), size = (1200, 600)
)
#jl savefig("tutorial_calculated_potential_undep.pdf") # hide
#md savefig("tutorial_calculated_potential_undep.pdf") # hide
#md savefig("tutorial_calculated_potential_undep.svg"); nothing # hide
#md # [![tutorial_calculated_potential_undep](tutorial_calculated_potential_undep.svg)](tutorial_calculated_potential_undep.pdf)

# 

is_depleted(sim_undep.point_types)

# Compare both volumes:

println("Depleted:   ", get_active_volume(sim.point_types))
println("Undepleted: ", get_active_volume(sim_undep.point_types));


# ## Electric field calculation

# Calculate the electric field of the fully depleted detector, given the already calculated electric potential:


calculate_electric_field!(sim, n_points_in_φ = 72)

plot(sim.electric_field, full_det = true, φ = 0.0, size = (700, 700))
plot_electric_fieldlines!(sim, full_det = true, φ = 0.0)
#jl savefig("tutorial_electric_field.pdf") # hide
#md savefig("tutorial_electric_field.pdf") # hide
#md savefig("tutorial_electric_field.svg"); nothing # hide
#md # [![tutorial_electric_field](tutorial_electric_field.svg)](tutorial_electric_field.pdf)


# ## Simulation of charge drifts

# Any charge drift model can be used for the calculation of the electric field. If no model is explicitely given, the `ElectricFieldChargeDriftModel` is used. Other configurations are saved in their configuration files and can be found under:

# `<package_directory>/examples/example_config_files/ADLChargeDriftModel/<config_filename>.yaml.`

# Set the charge drift model of the simulation:

charge_drift_model = ADLChargeDriftModel()
sim.detector = SolidStateDetector(sim.detector, charge_drift_model)

# Now, let's create an "random" multi-site event:

starting_positions = [ CylindricalPoint{T}( 0.020, deg2rad(10), 0.015 ),
                       CylindricalPoint{T}( 0.015, deg2rad(20), 0.045 ),
                       CylindricalPoint{T}( 0.022, deg2rad(35), 0.025 ) ]
energy_depos = T[1460, 609, 1000] * u"keV" # are needed later in the signal generation

evt = Event(starting_positions, energy_depos);

time_step = 5u"ns"
drift_charges!(evt, sim, Δt = time_step)

plot(sim.detector, size = (700, 700))
plot!(evt.drift_paths)
#jl savefig("tutorial_drift_paths.pdf") # hide
#md savefig("tutorial_drift_paths.pdf") # hide
#md savefig("tutorial_drift_paths.svg"); nothing # hide
#md # [![tutorial_drift_paths](tutorial_drift_paths.svg)](tutorial_drift_paths.pdf)


# ## Weighting potential calculation

# We need weighting potentials to simulate the detector charge signal induced by drifting charges. We'll calculate the weighting potential for the point contact and the outer shell of the detector:

for contact in sim.detector.contacts
    calculate_weighting_potential!(sim, contact.id, refinement_limits = [0.2, 0.1, 0.05, 0.01], n_points_in_φ = 2, verbose = false)
end

plot(
    plot(sim.weighting_potentials[1]),
    plot(sim.weighting_potentials[2]),
    size = (900, 700)
)
#jl savefig("tutorial_weighting_potentials.pdf") # hide
#md savefig("tutorial_weighting_potentials.pdf") # hide
#md savefig("tutorial_weighting_potentials.svg"); nothing # hide
#md # [![tutorial_weighting_potentials](tutorial_weighting_potentials.svg)](tutorial_weighting_potentials.pdf)

# ## Detector Capacitance Matrix

# After the weighting potentials are calculated, the detector capacitance matrix can be calculated in the **Maxwell Capacitance Matrix Notation**:

calculate_capacitance_matrix(sim)

# See [Capacitances](@ref) for more information.

# ## Detector waveform generation

# Given an interaction at an arbitrary point in the detector, we can now simulate the charge drift and the resulting detector charge signals (e.g. at the point contact):

simulate!(evt, sim) # drift_charges + signal generation of all channels

p_pc_signal = plot( evt.waveforms[1], lw = 1.5, xlims = (0, 1100), xlabel = "Time", unitformat = :slash,
                    legend = false, tickfontsize = 12, ylabel = "Charge", guidefontsize = 14)
#jl savefig("tutorial_waveforms.pdf") # hide
#md savefig("tutorial_waveforms.pdf") # hide
#md savefig("tutorial_waveforms.svg"); nothing # hide
#md # [![tutorial_waveforms](tutorial_waveforms.svg)](tutorial_waveforms.pdf)

# SolidStateDetectors.jl also allows to separate the waveform into the two contributions from electrons and holes

contact_id = 1
plot_electron_and_hole_contribution(evt, sim, contact_id, xlims = (0, 1100), xlabel = "Time",
                    legend = :topleft, tickfontsize = 12, ylabel = "Charge", guidefontsize = 14)
#jl savefig("tutorial_waveform_contributions.pdf") # hide
#md savefig("tutorial_waveform_contributions.pdf") # hide
#md savefig("tutorial_waveform_contributions.svg"); nothing # hide
#md # [![tutorial_waveform_contributions](tutorial_waveform_contributions.svg)](tutorial_waveform_contributions.pdf)