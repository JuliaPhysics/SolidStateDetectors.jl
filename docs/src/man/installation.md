# Installation

This package is a registered package.

Install via

```julia
using Pkg; pkg"add SolidStateDetectors"
```

## Visualization / Plotting (Optional)

This package provides serveral [plot recipes](https://docs.juliaplots.org/latest/recipes/) for different outputs for the plotting package [Plots.jl](https://github.com/JuliaPlots/Plots.jl/).

In order to use these also install the [Plots.jl](https://github.com/JuliaPlots/Plots.jl/) package via

```julia
using Pkg; pkg"add Plots"
```

Load the [Plots.jl](https://github.com/JuliaPlots/Plots.jl/) package (and optionally the backend `pyplot`) via

```julia
using Plots
```

The backends supported by SolidStateDetectors.jl are `gr` and `pyplot`.
By default, `gr` is loaded when importing `Plots`.

For more information about the plot recipes of this package look up the [Plotting](@ref) section.

This documentation was build with
```@example
using Pkg, Plots # hide
pkgversion(m::Module) = Pkg.TOML.parsefile(joinpath(dirname(string(first(methods(m.eval)).file)), "..", "Project.toml"))["version"] # hide
Plots_version = pkgversion(Plots) # hide
GR_version = pkgversion(GR) # hide
print("Plots: v$(Plots_version) - GR: v$(GR_version)") # hide
```

## GPU Support in Field Calculations

The [Electric Potential](@ref) and individual [Weighting Potentials](@ref) can also be calculated on GPUs. SolidStateDetectors.jl uses [KernelAbstractions.jl](https://github.com/JuliaGPU/KernelAbstractions.jl) for GPU support.

In order to use your GPU, the Julia Packages [CUDA.jl](https://github.com/JuliaGPU/CUDA.jl) or [AMDGPU.jl](https://github.com/JuliaGPU/AMDGPU.jl) have to be installed and loaded. SolidStateDetectors.jl may run on other GPU types supported by KernelAbstractions.jl, but has not been tested with [oneAPI.jl](https://github.com/JuliaGPU/oneAPI.jl) or [Metal.jl](https://github.com/JuliaGPU/Metal.jl) yet.

For NVIDIA GPUs use:

```julia
using CUDA, SolidStateDetectors
# *Only* for KernelAbstractions < v0.9:
# using CUDAKernels
```

For AMD GPUs use:

```julia
using AMDGPU, SolidStateDetectors
# *Only* for KernelAbstractions < v0.9:
# using ROCKernels
```

Then, in any field calculation ([`calculate_electric_potential!`](@ref), [`calculate_weighting_potential!`](@ref), [`simulate!(::Simulation)`](@ref)) the keyword `device_array_type` can be set to choose the device on which the calculations should be performed.
The possibilities are:
```julia
device_array_type = Array # -> CPU (default)
device_array_type = CuArray # -> NVIDIA GPU
device_array_type = ROCArray # -> AMD GPU
```

### Example (NVIDIA)

```julia
using CUDA, SolidStateDetectors
# *Only* for KernelAbstractions < v0.9:
# using CUDAKernels

sim = Simulation(SSD_examples[:CGD])
calculate_electric_potential!( 
    sim, 
    device_array_type = CuArray, 
    refinement_limits = [0.2, 0.1, 0.05],
    depletion_handling = true
)
```

!!! note
    The AMD backend was not yet tested due to lack of an AMD GPU (we are working on that).
