using SolidStateDetectors
using Geant4
# using Plots
using Unitful
using JLD2

include("snapshotting.jl")
T = Float32

function build_detsim(::Type{T}, example_geometry::Symbol) where T<:Real
    detector_config_filename = SSD_examples[example_geometry]
    sim = Simulation{T}(detector_config_filename)
    calculate_electric_potential!(sim, convergence_limit = 1e-6, refinement_limits = [0.2, 0.1, 0.05, 0.01])
    calculate_electric_field!(sim, n_points_in_φ = 72)
    for contact in sim.detector.contacts
        calculate_weighting_potential!(sim, contact.id, refinement_limits = [0.2, 0.1, 0.05, 0.01], n_points_in_φ = 2, verbose = false)
    end
    charge_drift_model = ADLChargeDriftModel()
    sim.detector = SolidStateDetector(sim.detector, charge_drift_model);
    return sim
end

@snapshot sim = build_detsim(Float32, :InvertedCoax)

function make_mcevents(sim)
    source_1 = MonoenergeticSource(
        "gamma",                              # Type of particle beam
        2.615u"MeV",                          # Energy of particle
        CartesianPoint(0.065, 0., 0.05),      # Location of the source
        CartesianVector(-1,0,0),              # Direction of the source
        10u"°"                                # Opening angle of the source emission
    )

    app = G4JLApplication(sim, source_1, verbose = false);
    N_events = 50000
    mcevents = run_geant4_simulation(app, N_events)
end

@snapshot mcevents = make_mcevents(sim)