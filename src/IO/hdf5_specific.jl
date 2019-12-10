function ssd_write(dest::AbstractString, s::Simulation)
    if isfile(dest) @warn "Destination `$dest` already exists. Overwriting..." end
    HDF5.h5open(dest, "w") do h5f
        LegendHDF5IO.writedata( h5f, "SSD_Simulation", NamedTuple(s)  )
    end       
end  


function ssd_read(src::AbstractString, ::Type{Simulation})
    HDF5.h5open(src, "r") do h5f
        Simulation(LegendHDF5IO.readdata(h5f, "SSD_Simulation"));
    end     
end  

