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

include("ConstructiveSolidGeometry/ConstructiveSolidGeometry.jl")
using .ConstructiveSolidGeometry
using .ConstructiveSolidGeometry:
            CylindricalPoint, CartesianPoint, AbstractCoordinatePoint, _convert_point,
            CartesianVector, CylindricalVector, AbstractCoordinateVector,
            Cartesian, Cylindrical, AbstractCoordinateSystem, CoordinateSystemType,
            CartesianTicksTuple, CylindricalTicksTuple,
            Geometry, AbstractGeometry, AbstractSurfacePrimitive,
            parse_rotation_matrix, parse_translate_vector, parse_CSG_transformation,
            transform, CSG_dict, Transformations, combine_transformations,
            ConfigFileError, _parse_value
        
import .ConstructiveSolidGeometry: sample
export CartesianPoint, CartesianVector, CylindricalPoint

import Clustering
import DataStructures
import Distributions
import Tables
import TypedTables

import Base: size, sizeof, length, getindex, setindex!, axes, getproperty, broadcast,
             range, ndims, eachindex, enumerate, iterate, IndexStyle, eltype, in, convert,
             show, print, println, display, +, -, &

export SolidStateDetector
export SSD_examples

export Grid

export ElectricPotential, PointTypes, EffectiveChargeDensity, DielectricDistribution, WeightingPotential, ElectricField
export apply_initial_state!
export calculate_electric_potential!, calculate_weighting_potential!, calculate_electric_field!, calculate_drift_fields!
export ElectricFieldChargeDriftModel, ADLChargeDriftModel
export get_active_volume
export simulate_waveforms
export Simulation, simulate!
export Event, drift_charges!

const SSDFloat = Union{Float16, Float32, Float64}


include("Units.jl")
include("Axes/DiscreteAxis.jl")
include("World/World.jl")
include("Grids/Grids.jl")

include("Types/Types.jl")

include("MaterialProperties/MaterialProperties.jl")
include("Config/Config.jl")
include("ChargeDensities/ChargeDensities.jl")
include("ImpurityDensities/ImpurityDensities.jl")
include("ChargeDriftModels/ChargeDriftModels.jl")
include("SolidStateDetector/DetectorGeometries.jl")

include("PotentialSimulation/PotentialSimulation.jl")

include("ElectricField/ElectricField.jl")

#include("ChargeCloudModels/ChargeCloudModels.jl")
include("ChargeDrift/ChargeDrift.jl")
include("SignalGeneration/SignalGeneration.jl")

include("ChargeStatistics/ChargeStatistics.jl")
include("ChargeClustering/ChargeClustering.jl")

include("Simulation/Simulation.jl")
include("Event/Event.jl")
include("MCEventsProcessing/MCEventsProcessing.jl")

include("IO/IO.jl")

include("examples.jl")

include("PlotRecipes/PlotRecipes.jl")

function __init__()
    @require HDF5="f67ccb44-e63f-5c2f-98bd-6dc0ccc4ba2f" begin
        @require LegendHDF5IO="c9265ca6-b027-5446-b1a4-febfa8dd10b0" begin
            include("IO/hdf5_specific.jl")
        end
        include("MCEventsProcessing/MCEventsProcessing_hdf5.jl")
    end
end

end # module
