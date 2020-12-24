# This file is a part of SolidStateDetectors.jl, licensed under the MIT License (MIT).

__precompile__(true)

module SolidStateDetectors

using LinearAlgebra
using Random
using Statistics

using ArraysOfArrays
using FillArrays
using Formatting
using Interpolations
using IntervalSets
using JSON
using LaTeXStrings
using ParallelProcessingTools
using ProgressMeter
using RadiationDetectorSignals
using RecipesBase
using Requires
using Rotations
using StaticArrays
using StatsBase
using Unitful
using YAML

import Clustering
import Distributions
import Tables
import TypedTables

import Base: size, sizeof, length, getindex, setindex!, axes, getproperty, broadcast,
             range, ndims, eachindex, enumerate, iterate, IndexStyle, eltype, in
import Base: show, print, println, display, +, -, &
import Base.convert

export SolidStateDetector
export SSD_examples

export Grid
# export CylindricalPoint, CartesianPoint

export ElectricPotential, PointTypes, EffectiveChargeDensity, DielectricDistribution, WeightingPotential, ElectricField
export apply_initial_state!
export calculate_electric_potential!, calculate_weighting_potential!, calculate_electric_field!
export update_till_convergence!, refine!
export set_charge_drift_model!, calculate_drift_fields!
export get_active_volume
export generate_charge_signals, generate_charge_signals!
export ElectricFieldChargeDriftModel, ADLChargeDriftModel
export Simulation, simulate!
export Event, drift_charges!

const SSDFloat = Union{Float16, Float32, Float64}

struct ConfigFileError <: Exception
    msg::AbstractString
end
Base.showerror(io::IO, e::ConfigFileError) = print(io, "ConfigFileError: ", e.msg)

include("ConstructiveSolidGeometry/ConstructiveSolidGeometry.jl")
# using .ConstructiveSolidGeometry
include("Geometries/Geometries.jl")

include("Axes/DiscreteAxis.jl")
include("World/World.jl")
include("Grids/Grids.jl")

include("Types/Types.jl")

include("MaterialProperties/MaterialProperties.jl")
include("Config/Config.jl")
include("ChargeDensityModels/ChargeDensityModels.jl")
include("SolidStateDetector/DetectorGeometries.jl")
include("GeometryRounding.jl")

include("PotentialSimulation/PotentialSimulation.jl")

include("ElectricField/ElectricField.jl")

include("ChargeDriftModels/ChargeDriftModels.jl")
include("ChargeCloudModels/ChargeCloudModels.jl")
include("ChargeDrift/ChargeDrift.jl")
include("SignalGeneration/SignalGeneration.jl")

include("ChargeStatistics/ChargeStatistics.jl")
include("ChargeClustering/ChargeClustering.jl")

include("Simulation/Simulation.jl")
include("Event/Event.jl")
include("MCEventsProcessing/MCEventsProcessing.jl")

include("IO/IO.jl")

include("examples.jl")

include("plotting/plotting.jl")

function __init__()
    @require HDF5="f67ccb44-e63f-5c2f-98bd-6dc0ccc4ba2f" begin
        @require LegendHDF5IO="c9265ca6-b027-5446-b1a4-febfa8dd10b0" begin
            include("IO/hdf5_specific.jl")
        end
        include("MCEventsProcessing/MCEventsProcessing_hdf5.jl")
    end
end

end # module
