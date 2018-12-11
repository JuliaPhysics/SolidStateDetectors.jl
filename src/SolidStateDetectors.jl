# This file is a part of SolidStateDetectors.jl, licensed under the MIT License (MIT).

__precompile__(true)

module SolidStateDetectors

using LinearAlgebra
using Random
using Statistics

using ArraysOfArrays
using CoordinateTransformations
using Interpolations
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

import Base: size, getindex, length
import Base: show, print, println, +

import Plots: grid, @layout

export SolidStateDetector
export SSD_examples

export calculate_electric_potential, calculate_weighting_potential

export AbstractChargeDriftModels, get_electron_drift_field, get_hole_drift_field
export VacuumChargeDriftModel, ADLChargeDriftModel
export get_active_volume
export ElectricPotential, CylindricalGrid, PointTypes
export generate_charge_signals!, generate_charge_signals


include("Types/Types.jl")
include("MaterialProperties/MaterialProperties.jl")
include("Geometries/Geometries.jl")
include("DetectorGeometries/DetectorGeometries.jl")
include("Grids/Grids.jl")
include("Potentials/Potentials.jl")
include("ElectricField/ElectricField.jl")
include("ChargeDriftModels/ChargeDriftModels.jl")
include("ChargeDrift/ChargeDrift.jl")
include("ChargeStatistics/ChargeStatistics.jl")
include("ChargeClustering/ChargeClustering.jl")

include("IO/IO.jl")

include("examples.jl")

end # module
