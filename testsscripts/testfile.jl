using SolidStateDetectors
using LegendHDF5IO

events_in = lh5open("testsscripts/simulation_output.lh5", "r") do h
    LegendHDF5IO.readdata(h.data_store, "SimulationData")
end

include("sim.jl")

starting_positionsv = [CartesianVector{T}(-0.02, 0.015, 0.04), 
    CartesianVector{T}(0.015, -0.012, 0.02), 
    CartesianVector{T}(0.01, -0.025, 0.01)]
energy_depos = T[1460, 609, 1000] * u"keV"

