function ssd_write(filename::AbstractString, sim::Simulation)
    if isfile(filename) @warn "Destination `$filename` already exists. Overwriting..." end
    HDF5.h5open(filename, "w") do h5f
        LegendHDF5IO.writedata( h5f, "SSD_Simulation", NamedTuple(sim)  )
    end       
end  

function ssd_read(filename::AbstractString, ::Type{Simulation})
    HDF5.h5open(filename, "r") do h5f
        Simulation(LegendHDF5IO.readdata(h5f, "SSD_Simulation"));
    end     
end  

