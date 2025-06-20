using Random
using SolidStateDetectors
using StructArrays, ArraysOfArrays
using Unitful
using RadiationDetectorSignals
using StructArrays, TypedTables
using JLD2

using SolidStateDetectors: internal_length_unit


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


function generate_charges(;
    n_events::Integer = 100, n_hits_avg::Integer = 10, n_charges_avg::Integer = 100,
    hit_dist_scale = T(50) * u"mm", charge_scale = T(1000) * u"keV", charge_cloud_scale = T(2) * u"mm"
)
    flat_evtno = Vector{Int}()
    flat_edep = Vector{typeof(charge_scale)}()
    PT = typeof(hit_dist_scale)
    flat_pos = StructVector{CartesianPoint{PT}}(undef, 0)

    for evtno in 1:n_events
        n_hits = rand(round(Int, 0.5 * n_hits_avg):round(Int, 1.5 * n_hits_avg))
        for _ in 1:n_hits
            n_charges = rand(round(Int, 0.5 * n_charges_avg):round(Int, 1.5 * n_charges_avg))

            pos_center = cartesian_zero + hit_dist_scale/2 * CartesianVector(2*rand(T)-1, 2*rand(T)-1, 2*rand(T)-1)

            edep = charge_scale * randexp(T, n_charges)

            pos = Ref(pos_center) .+ charge_cloud_scale .* CartesianVector.(randn(T, n_charges), randn(T, n_charges), randn(T, n_charges))

            idxs = let semi = sim.detector.semiconductor
                findall(in.(ustrip.(internal_length_unit, pos), Ref(semi)))
            end

            filtered_edep = edep[idxs]
            filtered_pos = pos[idxs]

            append!(flat_evtno, fill(evtno, length(filtered_edep)))
            append!(flat_edep, filtered_edep)
            append!(flat_pos, filtered_pos)
        end
    end
    return (;flat_evtno, flat_edep, flat_pos)
end

@snapshot generated_charges = generate_charges()

(;flat_evtno, flat_edep, flat_pos) = generated_charges


#=
using Plots
plot!(flat_pos[1:1000])
=#

flat_sa = StructArray(
    evtno = flat_evtno,
    edep = flat_edep,
    pos = flat_pos,
)

nested_sa = StructArray{typeof(flat_sa)}(
    evtno = consgroupedview(flat_evtno, flat_evtno),
    edep = consgroupedview(flat_evtno, flat_edep),
    pos = consgroupedview(flat_evtno, flat_pos),
)


mcevents = Table(
    evtno = consgroupedview(flat_evtno, flat_evtno),
    edep = consgroupedview(flat_evtno, flat_edep),
    pos = consgroupedview(flat_evtno, flat_pos),
)

# JLD2.save("../test_sim.jld2", "sim", sim)
# JLD2.save("../test_mcevents.jld2", "mcevents", mcevents)