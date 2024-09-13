include("SigGenInterface.jl")
include("ParseConfigFiles.jl")

NamedTuple(::Missing) = (object = "missing",)
NamedTuple(d::Dict) = (dict_json_string = json(d),)
Base.convert(T::Type{NamedTuple}, x::Dict) = T(x)
Dict(nt::NamedTuple) = JSON.parse(nt.dict_json_string)


"""
    ssd_write(filename::AbstractString, sim::Simulation)
    
Converts a [`Simulation`](@ref) to a `NamedTuple` and writes it to a HDF5 file 
with a given `filename` using [LegendHDF5IO.jl](https://github.com/legend-exp/LegendHDF5IO.jl).

## Arguments
* `filename::AbstractString`: Filename of the HDF5 file.
* `sim::Simulation`: [`Simulation`](@ref) that should be written to the HDF5 file.

## Example 
```julia
using LegendHDF5IO
using SolidStateDetectors
sim = Simulation(SSD_examples[:InvertedCoax])
simulate!(sim)
ssd_write("example_sim.h5", sim)
```

!!! warn
    If a file with `filename` already exists, it will be overwritten by this method.

!!! note 
    In order to use this method, the package [LegendHDF5IO.jl](https://github.com/legend-exp/LegendHDF5IO.jl) 
    has to be loaded before loading SolidStateDetectors.jl.

See also [`ssd_read`](@ref).
"""
function ssd_write end
export ssd_write


"""
    ssd_read(filename::AbstractString, ::Type{Simulation})
    
Reads a [`Simulation`](@ref) from a HDF5 file with a given `filename` 
using [LegendHDF5IO.jl](https://github.com/legend-exp/LegendHDF5IO.jl).

## Arguments
* `filename::AbstractString`: Filename of the HDF5 file.

## Example 
```julia
using LegendHDF5IO
using SolidStateDetectors
sim = ssd_read("example_sim.h5", Simulation)
```

!!! note 
    In order to use this method, the package [LegendHDF5IO.jl](https://github.com/legend-exp/LegendHDF5IO.jl) 
    has to be loaded before loading SolidStateDetectors.jl.

See also [`ssd_write`](@ref).
"""
function ssd_read end
export ssd_read


