# This file is a part of SolidStateDetectors.jl, licensed under the MIT License (MIT).

__precompile__(true)

module SolidStateDetectors

using LinearAlgebra
using Random
using Statistics

using ArraysOfArrays
using Interpolations
using IntervalSets
using JSON
using LaTeXStrings
using NamedTupleTools
using ParallelProcessingTools
using ProgressMeter
using RecipesBase
using StaticArrays
using Unitful
using YAML

import Clustering
import Distributions
import Tables
import TypedTables

import Base: size, sizeof, length, getindex, setindex!, axes, range, ndims, eachindex, enumerate, iterate, IndexStyle, eltype, in
import Base: show, print, println, display, +, -, &
import Base.convert

const SSD = SolidStateDetectors; export SSD
export SolidStateDetector
export SSD_examples

export Grid, CylindricalPoint, CartesianPoint

export ElectricPotential, PointTypes, ChargeDensity, DielectricDistribution, WeightingPotential, ElectricField
export calculate_electric_potential!, calculate_weighting_potential!, get_active_volume
export generate_charge_signals, generate_charge_signals!
export AbstractChargeDriftModel
export VacuumChargeDriftModel, ADLChargeDriftModel
export Simulation, simulate!

## temporary exports for easier debugging
export cyp, cap # CylindricalPoint, CartesianPoint
export point_type
export geom_round
export get_velocity_vector
export get_crossing_pos


const SSDFloat = Union{Float16, Float32, Float64}

struct ConfigFileError <: Exception 
    msg::AbstractString
end
Base.showerror(io::IO, e::ConfigFileError) = print(io, "ConfigFileError: ", e.msg)

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
include("ChargeDrift/ChargeDrift.jl")
include("SignalGeneration/SignalGeneration.jl")

include("ChargeStatistics/ChargeStatistics.jl")
include("ChargeClustering/ChargeClustering.jl")

include("Simulation/Simulation.jl")
include("Event/Event.jl")

include("IO/IO.jl")

include("examples.jl")

end # module
