# IO

After simulating the potentials and fields of a detector setup, the results should be saved to a file to avoid recalculating them every time the user starts a program.

## Saving output with JLD2

One easy way to do this is using [JLD2.jl](https://github.com/JuliaIO/JLD2.jl) and [FileIO.jl](https://github.com/JuliaIO/FileIO.jl).

Simulation results can be saved to a JLD file using `FileIO.save`:
```julia
using SolidStateDetectors
sim = Simulation("<config-file-name>")
# ...

using FileIO
FileIO.save("<name-of-simulation-file>.jld", Dict("Simulation" => sim))
```

It can be read back in using `FileIO.load`:
```julia
using FileIO
sim = FileIO.load("<name-of-simulation.file>.jld", "Simulation")
```

## Saving output with HDF5

One more compact way of saving simulation results is based on converting the output to a `NamedTuple` and saving it to a HDF5 file. This requires the (unregistered) service packages [LegendDataTypes.jl](https://github.com/legend-exp/LegendHDF5IO.jl) and [LegendHDF5IO.jl](https://github.com/legend-exp/LegendHDF5IO.jl).

Install the required packages once by running:
```julia
import Pkg
Pkg.add(url="https://github.com/legend-exp/LegendDataTypes.jl.git")
Pkg.add(url="https://github.com/legend-exp/LegendHDF5IO.jl.git")
```

Simulation output can be written to a HDF5 file using [`ssd_write`](@ref):
```julia
using SolidStateDetectors 
using LegendHDF5IO
# ...

ssd_write("<name-of-simulation-file>.h5", sim)
```

The data stored in the HDF5 file can be read using [`ssd_read`](@ref):
```julia
using SolidStateDetectors
using LegendHDF5IO
ssd_read("<name-of-simulation-file>.h5", Simulation)
```