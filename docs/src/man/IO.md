# IO

Right now (electric & weighting) potentials and point types can easily be saved and loaded via [HDF5.jl](https://github.com/JuliaIO/HDF5.jl).

```julia
using SolidStateDetectors, HDF5
detector = SolidStateDetector(SSD_examples[:InvertedCoax])
E_pot, point_types = calculate_electric_potential(detector)

# write to HDF5
h5f = h5open("InvertedCoaxSimulation.hdf5", "w") 
g_E_pot = g_create(h5f, "Electric Potential")
g_point_types = g_create(h5f, "Point Types")
SolidStateDetectors.write_to_hdf5(g_E_pot, E_pot) 
SolidStateDetectors.write_to_hdf5(g_point_types, point_types) 
close(h5f)
```

```julia
using SolidStateDetectors, HDF5

# read from HDF5
h5f = h5open("InvertedCoaxSimulation.hdf5", "r") 
g_E_pot = g_open(h5f, "Electric Potential")
g_point_types = g_open(h5f, "Point Types")
E_pot = SolidStateDetectors.read_from_hdf5(g_E_pot, ElectricPotential)
point_types = SolidStateDetectors.read_from_hdf5(g_point_types, PointTypes)
close(h5f)
```
