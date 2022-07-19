# This file is a part of SolidStateDetectors.jl, licensed under the MIT License (MIT).

__precompile__(true)

module SolidStateDetectors

using LinearAlgebra
using Random
using Statistics

using Adapt
using ArraysOfArrays
using FillArrays
using Formatting
using GPUArrays
using Interpolations
using IntervalSets
using JSON
using KernelAbstractions
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
using UnitfulAtomic
using UnitfulRecipes
using YAML

include("ConstructiveSolidGeometry/ConstructiveSolidGeometry.jl")
using .ConstructiveSolidGeometry
using .ConstructiveSolidGeometry:
            CylindricalPoint, CartesianPoint, AbstractCoordinatePoint, _convert_point,
            CartesianVector, CylindricalVector, Plane, AbstractCoordinateVector,
            Cartesian, Cylindrical, AbstractCoordinateSystem, CoordinateSystemType,
            Geometry, AbstractGeometry, AbstractSurfacePrimitive,
            parse_rotation_matrix, parse_translate_vector, parse_CSG_transformation,
            transform, CSG_dict, Transformations, combine_transformations,
            ConfigFileError, _parse_value,
            LengthQuantity, AngleQuantity, get_scale
        
import .ConstructiveSolidGeometry: sample, to_internal_units, from_internal_units
export CartesianPoint, CartesianVector, CylindricalPoint

import Clustering
import DataStructures
import Distributions
import IntervalSets
import Tables
import TypedTables

import Base: size, sizeof, length, getindex, setindex!, axes, getproperty, broadcast,
             range, ndims, eachindex, enumerate, iterate, IndexStyle, eltype, in, convert,
             show, print, println, display, +, -, &

export SolidStateDetector
export SSD_examples

export Grid, GridPoint

export ElectricPotential, PointTypes, EffectiveChargeDensity, DielectricDistribution, WeightingPotential, ElectricField
export apply_initial_state!
export calculate_electric_potential!, calculate_weighting_potential!, calculate_electric_field!, calculate_drift_fields!
export ElectricFieldChargeDriftModel, ADLChargeDriftModel
export get_active_volume, is_depleted
export calculate_stored_energy, calculate_mutual_capacitance, calculate_capacitance_matrix
export simulate_waveforms
export Simulation, simulate!
export Event, drift_charges!
export add_baseline_and_extend_tail
export NBodyChargeCloud

const SSDFloat = Union{Float16, Float32, Float64}


include("Units.jl")
include("Axes/DiscreteAxis.jl")
include("World/World.jl")
include("Symmetry/Symmetry.jl")
include("Grids/Grids.jl")

include("MaterialProperties/MaterialProperties.jl")
include("Config/Config.jl")
include("ChargeDensities/ChargeDensities.jl")
include("ImpurityDensities/ImpurityDensities.jl")
include("ChargeDriftModels/ChargeDriftModels.jl")
include("SolidStateDetector/DetectorGeometries.jl")

include("ScalarPotentials/ScalarPotential.jl")
include("PotentialCalculation/PotentialCalculation.jl")

include("ElectricField/ElectricField.jl")

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

include("PlotRecipes/PlotRecipes.jl")
export @P_str # protected strings to overwrite plot labels with units

function __init__()
    @require HDF5="f67ccb44-e63f-5c2f-98bd-6dc0ccc4ba2f" begin
        @require LegendHDF5IO="c9265ca6-b027-5446-b1a4-febfa8dd10b0" begin
            include("IO/hdf5_specific.jl")
        end
        include("MCEventsProcessing/MCEventsProcessing_hdf5.jl")
    end
end

end # module
