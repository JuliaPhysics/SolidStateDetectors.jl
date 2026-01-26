# This file is a part of SolidStateDetectors.jl, licensed under the MIT License (MIT).

module SolidStateDetectorsLegendHDF5IOExt

import ..LegendHDF5IO

using ..SolidStateDetectors
using ..SolidStateDetectors: RealQuantity, SSDFloat, to_internal_units, chunked_ranges, LengthQuantity
using ..SolidStateDetectors.ConstructiveSolidGeometry: AbstractCoordinatePoint

using ArraysOfArrays
import Tables
using TypedTables, Unitful, UnitfulAtomic
using Format

function SolidStateDetectors.ssd_write(filename::AbstractString, sim::Simulation)
    if isfile(filename) @warn "Destination `$filename` already exists. Overwriting..." end
    LegendHDF5IO.lh5open(filename, "w") do h5f
        LegendHDF5IO.writedata(h5f.data_store, "SSD_Simulation", NamedTuple(sim)  )
    end       
end  

function SolidStateDetectors.ssd_read(filename::AbstractString, ::Type{Simulation})
    LegendHDF5IO.lh5open(filename, "r") do h5f
        Simulation(LegendHDF5IO.readdata(h5f.data_store, "SSD_Simulation"));
    end     
end  



function SolidStateDetectors.simulate_waveforms( mcevents::AbstractVector{<:NamedTuple}, sim::Simulation{T},
                             output_dir::AbstractString, 
                             output_base_name::AbstractString = "generated_waveforms";
                             chunk_n_physics_events::Int = 1000, 
                             Δt::RealQuantity = 4u"ns",
                             max_nsteps::Int = 1000,
                             diffusion::Bool = false,
                             self_repulsion::Bool = false,
                             number_of_carriers::Int = 1,
                             number_of_shells::Int = 1,
                             signal_unit::Unitful.Units = u"e_au",
                             max_interaction_distance::Union{<:Real, <:LengthQuantity} = NaN,
                             end_drift_when_no_field::Bool = true,
                             geometry_check::Bool = false,
                             verbose = false) where {T <: SSDFloat}
    n_total_physics_events = length(mcevents)
    Δtime = T(to_internal_units(Δt)) 
    n_contacts = length(sim.detector.contacts)
    @info "Detector has $(n_contacts) contact(s)"
    if !ispath(output_dir) mkpath(output_dir) end
    nfmt(i::Int) = format(i, zeropadding = true, width = length(digits(n_total_physics_events)))
    evt_ranges = chunked_ranges(n_total_physics_events, chunk_n_physics_events)
    @info "-> $(length(mcevents)) events to simulate."

    for evtrange in evt_ranges
        ofn = joinpath(output_dir, "$(output_base_name)_evts_$(nfmt(first(evtrange)))-$(nfmt(last(evtrange))).h5")
        @info "Now simulating $(evtrange) and storing it in\n\t \"$ofn\""
        mcevents_sub = simulate_waveforms(mcevents[evtrange], sim; Δt, max_nsteps, diffusion, self_repulsion, number_of_carriers, number_of_shells, signal_unit, max_interaction_distance, end_drift_when_no_field, geometry_check, verbose)

        # # LH5 couldn't handle CartesianPoint, turn positions into CartesianVectors which will be saved as SVectors
        # pos_vec = VectorOfVectors([[CartesianPoint(p...) - cartesian_zero for p in ps] for ps in mcevents_sub.pos])
        # new_mcevents_sub = TypedTables.Table(merge(Tables.columns(mcevents_sub), (pos = pos_vec,)))

        LegendHDF5IO.lh5open(ofn, "w") do h5f
            LegendHDF5IO.writedata(h5f.data_store, "generated_waveforms", mcevents_sub)
        end        
    end    
end


# Add support for reading and writing CartesianPoint and CylindricalPoint using LegendHDF5IO
function __init__()
    LegendHDF5IO._datatype_dict["cartesianpoint"] = CartesianPoint
    LegendHDF5IO._datatype_dict["cylindricalpoint"] = CylindricalPoint
end

LegendHDF5IO.datatype_to_string(::Type{<:CartesianPoint}) = "cartesianpoint"
LegendHDF5IO.datatype_to_string(::Type{<:CylindricalPoint}) = "cylindricalpoint"



# support writing points to HDF5 file using LH5Array
function LegendHDF5IO.create_entry(parent::LegendHDF5IO.LHDataStore, name::AbstractString, 
    pt::T; kwargs...) where {T <: AbstractCoordinatePoint}
    LegendHDF5IO.create_entry(parent, name, getindex.(Ref(pt), 1:3); kwargs...)
    LegendHDF5IO.setdatatype!(parent.data_store[name], T)
    return nothing
end

function LegendHDF5IO.LH5Array(ds::LegendHDF5IO.HDF5.Dataset, T::Type{<:AbstractCoordinatePoint})
    T(read(ds) .* LegendHDF5IO.getunits(ds)...)
end

# support writing arrays of points to HDF5 file using LH5Array
function LegendHDF5IO.create_entry(parent::LegendHDF5IO.LHDataStore, name::AbstractString, 
    pts::T; kwargs...) where {T <: AbstractVector{<:AbstractCoordinatePoint}}
    LegendHDF5IO.create_entry(parent, name, VectorOfVectors([getindex.(Ref(pt), 1:3) for pt in pts]); kwargs...)
    LegendHDF5IO.setdatatype!(parent.data_store[name], T)
    return nothing
end

function LegendHDF5IO.LH5Array(ds::LegendHDF5IO.HDF5.Group, T::Type{<:AbstractVector{<:CartesianPoint}})
    data = LegendHDF5IO.LH5Array(ds, VectorOfVectors)
    Vector([CartesianPoint(d...) for d in data])
end

function LegendHDF5IO.LH5Array(ds::LegendHDF5IO.HDF5.Group, T::Type{<:AbstractVector{<:CylindricalPoint}})
    data = LegendHDF5IO.LH5Array(ds, VectorOfVectors)
    Vector([CylindricalPoint(d...) for d in data])
end

# support writing points to HDF5 file using readdata/writedata
function LegendHDF5IO.writedata(
    output::LegendHDF5IO.HDF5.H5DataStore, name::AbstractString,
    pt::AbstractCoordinatePoint,
    fulldatatype::DataType = typeof(pt)
)
    data = getindex.(Ref(pt), 1:3)
    LegendHDF5IO.writedata(output, name, data, fulldatatype)
end

function LegendHDF5IO.readdata(
    input::LegendHDF5IO.HDF5.H5DataStore, name::AbstractString,
    fulldatatype::Type{T}
) where {T <: AbstractCoordinatePoint}
    data = LegendHDF5IO.readdata(input, name, Vector)
    T(data...)
end


# support writing arrays of points to HDF5 file using readdata/writedata
function LegendHDF5IO.writedata(
    output::LegendHDF5IO.HDF5.H5DataStore, name::AbstractString,
    pts::AbstractVector{<:AbstractCoordinatePoint},
    fulldatatype::DataType = typeof(pts)
)
    data = VectorOfVectors([getindex.(Ref(pt), 1:3) for pt in pts])
    LegendHDF5IO.writedata(output, name, data, fulldatatype)
    LegendHDF5IO.setdatatype!(output[name], fulldatatype)
    return nothing
end

function LegendHDF5IO.readdata(
    input::LegendHDF5IO.HDF5.H5DataStore, name::AbstractString,
    fulldatatype::Type{T}
) where {T <:AbstractVector{<:CartesianPoint}}
    data = LegendHDF5IO.readdata(input, name, VectorOfVectors)
    output = [CartesianPoint(d...) for d in data]
    @assert output isa T
    return output
end

function LegendHDF5IO.readdata(
    input::LegendHDF5IO.HDF5.H5DataStore, name::AbstractString,
    fulldatatype::Type{T}
) where {T <:AbstractVector{<:CylindricalPoint}}
    data = LegendHDF5IO.readdata(input, name, VectorOfVectors)
    output = [CylindricalPoint(d...) for d in data]
    @assert output isa T
    return output
end

end # module LegendHDF5IO
