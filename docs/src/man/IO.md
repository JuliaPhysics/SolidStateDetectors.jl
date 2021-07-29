# IO

After simulating the potentials and fields of a detector setup, the results should be saved to a file to avoid recalculating them every time the user starts a program.

One easy way to do this is using [JLD2.jl](https://github.com/JuliaIO/JLD2.jl) and [FileIO.jl](https://github.com/JuliaIO/FileIO.jl).

Simulation results can be saved to a JLD file using
```yaml
using SolidStateDetectors
sim = Simulation("<config-file-name>")
# ...

using FileIO
FileIO.save("<name-of-simulation-file>.jld", Dict("Simulation" => sim)
```

It can be read back in using
```yaml 
using FileIO
sim = FileIO.load("<name-of-simulation.file>.jld", "Simulation")
```

Other more compact ways of saving simulation results are based on the HDF5 saving format and the package [HDF5.jl](https://github.com/JuliaIO/HDF5.jl).