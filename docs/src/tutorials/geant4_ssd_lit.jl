# # Geant4 Support

# SolidStateDetectors.jl provides an extension for Geant4.jl.
# The extension allows users to generate realistic event distributions resulting from particles emitted by a specified source, which could then in turn also be used as input for a waveform simulation.
#
# To use the extension, both `SolidStateDetectors` and `Geant4` have to be loaded.

using SolidStateDetectors
using Geant4

# In order to run Geant4 simulations, a `Geant4.G4JLApplication` needs to be defined, which needs to include both the detector geometry and the particle source.
# The extension features a function that creates a `Geant4.G4JLApplcation` object from a SSD `Simulation` object and a particle source.

using Plots
using Unitful

# Two types of particle source are pre-defined in `SolidStateDetectors`:
# 
# #### 1. `MonoenergeticSource`
# 
# This source emits particles with a fixed type and energy.

source_1 = MonoenergeticSource(
    "gamma",                              # Type of particle beam
    2.615u"MeV",                          # Energy of particle
    CartesianPoint(0.065, 0., 0.05),      # Location of the source
    CartesianVector(-1,0,0),              # Direction of the source
    10u"°"                                # Opening angle of the source emission
)

# - The particle type is given as a string (e.g. `"e-"` or `"gamma"`) and directly passed to `Geant4`. See the `Geant4` documentation on how to name the desired particle type.
# - The energy of the emitted particles is passed as a number with the required unit.
# - The `position` of the particle source relative to the origin is defined by a `CartesianPoint` (in units of `m`).
# - The source emits particles in a given `direction`, specified by a `CartesianVector` object. If the user provides no direction, the emission will be isotropic.
# - If an `opening_angle` is provided, the particles are emitted in a directed cone with the specified opening angle.

# #### 2. `IsotopeSource`
# 
# This source emits particles based on the radioactive decay chain of a given isotope.

source_2 = IsotopeSource(
    82,                                # Number of protons
    212,                               # Total number of nucleons
    0.0,                               # Excitation energy
    0.0,                               # Ion charge
    CartesianPoint(0.06, 0, 0.05),     # Location of the source
    CartesianVector(-1,0,0),           # Direction of the source
    10u"°"                             # Opening angle of the source emission
)

# The source is defined using
# - The number of protons `Z` and the number of nucleons `A` in the isotope.
# - The excitation energy
# - The charge of the isotope
# - The position, direction and opening angle from the source can be defined in the same way as for a `MonoenergeticSource`
# 
# The particle source can now be plotted together with the detector, as well as the direction in which particles are emitted.

T = Float32
sim = Simulation{T}(SSD_examples[:InvertedCoaxInCryostat])
plot(sim.detector, size = (500,500))
plot!(source_1)
#jl savefig("detector.pdf") # hide
#md savefig("detector.pdf") # hide
#md savefig("detector.svg"); nothing # hide
#md # [![detector](detector.svg)](detector.pdf) 


# A `Geant4.G4JLApplication` is built from a SSD `Simulation` `sim` and one of the previously defined particle sources, e.g. `source_1`.
# 
# Internally, the translated geometry description is written to a temporary GDML file. This file will in turn be read in by `Geant4.jl` and deleted afterwards. However, if needed, the resulting GDML file can also be saved by using the `Geant4.G4JLDetector(sim, "output_filename.gdml", save_gdml = true)` function.

app = G4JLApplication(sim, source_1, verbose = false);


# Having defined a `G4JLApplication`, it is now possible to generate events using the `run_geant4_simulation` function. The function will continue running until the desired number of energy depositions have been recorded.

N_events = 50000
events = run_geant4_simulation(app, N_events)

# Each entry of the table corresponds to one event and consists of five fields:
# - `evtno`: Event number
# - `detno`: Index of the detector where the energy was deposited
# - `thit`: Time of the interaction
# - `edep`: Amount of energy that was deposited in the detector
# - `pos`: Position where the interaction happened

# By extracting the position of each energy deposition from `events`, the spatial distribution of the events inside the detector can be visualized:

plot(sim.detector, show_passives = false, size = (500,500), fmt = :png)
plot!(source_1)
plot!(CartesianPoint.(broadcast(p -> ustrip.(u"m", p), events[1:1000].pos.data)), ms = 0.5, msw = 0, color=:black, label = "")
#jl savefig("events.pdf") # hide
#md savefig("events.pdf") # hide
#md savefig("events.svg"); nothing # hide
#md # [![events](events.svg)](events.pdf) 

# The output of `run_geant4_simulation` can be stored using the `LegendHDF5IO` package:

# ```Julia
# using LegendHDF5IO
# 
# lh5open("simulation_output.lh5", "w") do h
#     LegendHDF5IO.writedata(h.data_store, "SimulationData", events)
# end

# events_in = lh5open("simulation_output.lh5", "r") do h
#     LegendHDF5IO.readdata(h.data_store, "SimulationData")
# end
# ```

# In order to visualize the energy spectrum of the events in a histogram, the following code can be used:

using StatsBase

h = fit(Histogram, ustrip.(u"keV", sum.(events.edep)), Weights(fill(10,length(events.edep))), 0:10:3000)

plot(h, title = "Energy spectrum", bins = 500, yscale = :log10, st = :step, label = "")
xlims!(0,3000, xlabel = "E in keV", ylabel = "counts")
#jl savefig("spectrum.pdf") # hide
#md savefig("spectrum.pdf") # hide
#md savefig("spectrum.svg"); nothing # hide
#md # [![spectrum](spectrum.svg)](spectrum.pdf) 


# Now that a realistic event distribution has been generated, it can be used as input for the waveform simulation.
# To perform the waveform simulation, the detectors electric field and weighting potential have to be calculated first.

sim.detector = SolidStateDetector(sim.detector, ADLChargeDriftModel(T=T))
calculate_electric_potential!(sim, refinement_limits = [0.4,0.2,0.1,0.06], verbose = false)
calculate_electric_field!(sim, n_points_in_φ = 10)
calculate_weighting_potential!(sim, 1, refinement_limits = [0.4,0.2,0.1,0.06], verbose = false)


# The previously generated events from the Geant4 simulation can now be used as input for the `simulate_waveforms` function:

wf = simulate_waveforms(events[1:100], sim, Δt = 1u"ns", max_nsteps = 2000)
plot(wf[1:20].waveform, label = "")
#jl savefig("waveforms.pdf") # hide
#md savefig("waveforms.pdf") # hide
#md savefig("waveforms.svg"); nothing # hide
#md # [![waveforms](waveforms.svg)](waveforms.pdf) 



# We can add some baseline and tail to the pulses to match their lengths (in this case to 2000ns):

w = add_baseline_and_extend_tail.(wf.waveform, 100, 2000)
plot(w[1:20], label = "")
#jl savefig("wf_and_amplitude.pdf") # hide
#md savefig("wf_and_amplitude.pdf") # hide
#md savefig("wf_and_amplitude.svg"); nothing # hide
#md # [![wf_and_amplitude](wf_and_amplitude.svg)](wf_and_amplitude.pdf) 