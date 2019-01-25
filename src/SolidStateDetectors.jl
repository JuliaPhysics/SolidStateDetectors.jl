# This file is a part of SolidStateDetectors.jl, licensed under the MIT License (MIT).

__precompile__(true)

module SolidStateDetectors

using LinearAlgebra
using Random
using Statistics

using ArraysOfArrays
using CoordinateTransformations
using Interpolations
using IntervalSets
using JSON
using LaTeXStrings
using ParallelProcessingTools
using ProgressMeter
using RecipesBase
using StaticArrays
using Unitful

import Clustering
import Distributions
import Tables
import TypedTables

import Base: size, sizeof, length, getindex, setindex!, axes, range, ndims, eachindex, enumerate, iterate, IndexStyle, eltype
import Base: show, print, println, display, +

export SolidStateDetector
export SSD_examples

const SSD = SolidStateDetectors
export SSD

export calculate_electric_potential, calculate_weighting_potential

export AbstractChargeDriftModels, get_electron_drift_field, get_hole_drift_field
export VacuumChargeDriftModel, ADLChargeDriftModel
export get_active_volume
export Grid
export ElectricPotential, PointTypes, ChargeDensity, DielectricDistribution, WeightingPotential
export generate_charge_signals!, generate_charge_signals    

include("GeometryRounding.jl")

include("Axes/DiscreteAxis.jl")
include("Grids/Grids.jl")
include("Types/Types.jl")

include("MaterialProperties/MaterialProperties.jl")
include("Geometries/Geometries.jl")
include("DetectorGeometries/DetectorGeometries.jl")

include("Config/Config.jl")
include("PotentialSimulation/PotentialSimulation.jl")

include("ElectricField/ElectricField.jl")
include("ChargeDriftModels/ChargeDriftModels.jl")
include("ChargeDrift/ChargeDrift.jl")
include("ChargeStatistics/ChargeStatistics.jl")
include("ChargeClustering/ChargeClustering.jl")

include("IO/IO.jl")

include("examples.jl")

end # module
