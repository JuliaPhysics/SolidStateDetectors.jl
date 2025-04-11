# This file is a part of SolidStateDetectors.jl, licensed under the MIT License (MIT).

module SolidStateDetectorsLegendHDF5IOExt

import ..LegendHDF5IO

using SolidStateDetectors
using SolidStateDetectors: RealQuantity, SSDFloat, to_internal_units, chunked_ranges
using TypedTables, Unitful
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



function SolidStateDetectors.simulate_waveforms( mcevents::TypedTables.Table, sim::Simulation{T},
                             output_dir::AbstractString, 
                             output_base_name::AbstractString = "generated_waveforms";
                             chunk_n_physics_events::Int = 1000, 
                             Δt::RealQuantity = 4u"ns",
                             max_nsteps::Int = 1000,
                             diffusion::Bool = false,
                             self_repulsion::Bool = false,
                             number_of_carriers::Int = 1,
                             number_of_shells::Int = 1,
                             verbose = false) where {T <: SSDFloat}
    n_total_physics_events = length(mcevents)
    Δtime = T(to_internal_units(Δt)) 
    n_contacts = length(sim.detector.contacts)
    @info "Detector has $(n_contacts) contact(s)"
    contacts = sim.detector.contacts
    if !ispath(output_dir) mkpath(output_dir) end
    nfmt(i::Int) = format(i, zeropadding = true, width = length(digits(n_total_physics_events)))
    evt_ranges = chunked_ranges(n_total_physics_events, chunk_n_physics_events)
    @info "-> $(sum(length.(mcevents.edep))) energy depositions to simulate."

    for evtrange in evt_ranges
        ofn = joinpath(output_dir, "$(output_base_name)_evts_$(nfmt(first(evtrange)))-$(nfmt(last(evtrange))).h5")
        @info "Now simulating $(evtrange) and storing it in\n\t \"$ofn\""
        mcevents_sub = simulate_waveforms(mcevents[evtrange], sim; Δt, max_nsteps, diffusion, self_repulsion, number_of_carriers, number_of_shells, verbose)
      
        LegendHDF5IO.lh5open(ofn, "w") do h5f
            LegendHDF5IO.writedata(h5f.data_store, "generated_waveforms", mcevents_sub)
        end        
    end    
end

end # module LegendHDF5IO
