using SolidStateDetectors
using Geant4

using Plots
using Unitful

source_1 = MonoenergeticSource(
    "gamma",                              # Type of particle beam
    2.615u"MeV",                          # Energy of particle
    CartesianPoint(0.065, 0., 0.05),      # Location of the source
    CartesianVector(-1,0,0),              # Direction of the source
    10u"°"                                # Opening angle of the source emission
)

source_2 = IsotopeSource(
    82,                                # Number of protons
    212,                               # Total number of nucleons
    0.0,                               # Excitation energy
    0.0,                               # Ion charge
    CartesianPoint(0.06, 0, 0.05),     # Location of the source
    CartesianVector(-1,0,0),           # Direction of the source
    10u"°"                             # Opening angle of the source emission
)

T = Float32
sim = Simulation{T}(SSD_examples[:InvertedCoaxInCryostat])
plot(sim.detector, size = (500,500))
plot!(source_1)
savefig("detector.pdf") # hide

app = G4JLApplication(sim, source_1, verbose = false);

N_events = 50000
events = run_geant4_simulation(app, N_events)

plot(sim.detector, show_passives = false, size = (500,500), fmt = :png)
plot!(source_1)
plot!(CartesianPoint.(broadcast(p -> ustrip.(u"m", p), events[1:1000].pos.data)), ms = 0.5, msw = 0, color=:black, label = "")
savefig("events.pdf") # hide

using StatsBase

h = fit(Histogram, ustrip.(u"keV", sum.(events.edep)), Weights(fill(10,length(events.edep))), 0:10:3000)

plot(h, title = "Energy spectrum", bins = 500, yscale = :log10, st = :step, label = "")
xlims!(0,3000, xlabel = "E in keV", ylabel = "counts")
savefig("spectrum.pdf") # hide

sim.detector = SolidStateDetector(sim.detector, ADLChargeDriftModel(T=T))
calculate_electric_potential!(sim, refinement_limits = [0.4,0.2,0.1,0.06], verbose = false)
calculate_electric_field!(sim, n_points_in_φ = 10)
calculate_weighting_potential!(sim, 1, refinement_limits = [0.4,0.2,0.1,0.06], verbose = false)

wf = simulate_waveforms(events[1:100], sim, Δt = 1u"ns", max_nsteps = 2000)
plot(wf[1:20].waveform, label = "")
savefig("waveforms.pdf") # hide

w = add_baseline_and_extend_tail.(wf.waveform, 100, 2000)
plot(w[1:20], label = "")
savefig("wf_and_amplitude.pdf") # hide
