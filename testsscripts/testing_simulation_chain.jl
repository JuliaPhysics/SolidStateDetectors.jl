# include("gen_test_charges.jl")
using SolidStateDetectors
using SolidStateDetectors: drift_charges, _drift_charge!, Event, EHDriftPath, to_internal_units, Event, to_internal_units, interpolated_vectorfield, Electron, internal_length_unit
using Unitful
using BenchmarkTools
using JLD2

using ArraysOfArrays

SSD = SolidStateDetectors

T = Float32
sim = JLD2.load("../test_sim.jld2", "sim")

# ...
electric_field = interpolated_vectorfield(sim.electric_field)
time_step = 5u"ns"
max_nsteps = 1000
diffusion = false
self_repulsion = true
verbose = true

starting_positions = [CartesianPoint{T}(-0.02, 0.015, 0.04), 
    CartesianPoint{T}(0.015, -0.012, 0.02), 
    CartesianPoint{T}(0.01, -0.025, 0.01)]
energy_depos = T[1460, 609, 1000] * u"keV" # are needed later in the signal generation


evt2 = Event(starting_positions, energy_depos, 4)

time_step = 5u"ns"

start_pos2 = flatview(evt2.locations)
edep2 =  flatview(evt2.energies)
drift_paths::Vector{EHDriftPath{T}} = Vector{EHDriftPath{T}}(undef, length(start_pos2))
dt::T = T(to_internal_units(time_step))

drift_path_counter::Int = 0

n_hits::Int = length(start_pos2)
charges::Vector{T} = to_internal_units.(edep2) ./ to_internal_units(sim.detector.semiconductor.material.E_ionisation)

drift_path_e::Array{CartesianPoint{T},2} = Array{CartesianPoint{T},2}(undef, n_hits, max_nsteps)
timestamps_e::Vector{T} = Vector{T}(undef, max_nsteps)

@benchmark _drift_charge!(drift_path_e, timestamps_e, sim.detector, sim.point_types, sim.point_types.grid, start_pos2, -charges, dt, electric_field, Electron, diffusion=false, self_repulsion=true, verbose=true)
@benchmark for i in 1:20; _drift_charge!(drift_path_e, timestamps_e, sim.detector, sim.point_types, sim.point_types.grid, start_pos2, -charges, dt, electric_field, Electron, diffusion=false, self_repulsion=true, verbose=true); end

@profview for i in 1:20; _drift_charge!(drift_path_e, timestamps_e, sim.detector, sim.point_types, sim.point_types.grid, start_pos2, -charges, dt, electric_field, Electron, diffusion=false, self_repulsion=true, verbose=true); end





