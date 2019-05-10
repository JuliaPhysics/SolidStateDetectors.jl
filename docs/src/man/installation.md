# Installation

This package is a registered package.

Install via

```julia
using Pkg; pkg"add SolidStateDetectors"
```

## Vizualization / Plotting (Optional)

This package provides serveral [plot recipes](https://docs.juliaplots.org/latest/recipes/) for different outputs for the plotting package [Plots.jl](https://github.com/JuliaPlots/Plots.jl/).

In order to use these also install the [Plots.jl](https://github.com/JuliaPlots/Plots.jl/) package via

```julia
using Pkg; pkg"add Plots"
```

It is recommended to use `pyplot` as backend. Install via
```julia
using Pkg; pkg"add PyPlot"
```

Then you can load it via
```julia
using Plots; pyplot();
```