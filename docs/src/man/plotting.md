# Plotting

````@example tutorial
using Plots
using SolidStateDetectors
using Unitful

T = Float32
sim = Simulation{T}(SSD_examples[:BEGe]);
simulate!(sim, convergence_limit = 1e-6, refinement_limits = [0.2, 0.1, 0.05, 0.01]);
````

Besides the usual `Plots.jl`-keywords settings, like `size = (900, 900)`,
there are some additional settings to tune the plots which are described in the following sections.

## Detector Plots

The geometry of a detector, together with its environment, can be simply plotted via
````@example tutorial
det = sim.detector
plot(det, size = (900,900))
````

### Optional keywords:
* `seriestype`: Can be `:mesh3d` (default).
* `linewidth`: Sets the linewidth of the edges of the mesh.
* `linecolor`: Sets the linecolor of all edges of the mesh.
* `fillcolor`: Sets the face color of all faces.
* `fillalpha`: Sets the alpha value of all faces of the mesh.
* `show_semiconductor`: Whether also the primitives of the semiconductor should be plotted. Default is `false`
* `show_passives` = Whether also the primitives of the sourrinding objects of the detector should be plotted. Default is `true`

#### How does the plot recipe work?

The detector consists out of a semiconductor, `det.semiconductor`, its contacts, `sim.contacts`,
and its surrounding objects, `sim.passives`.
By default, see optional keywords, only the contacts and the surrounding objects are plotted
and the color of the contact primitives is defined internally through its `id`.

If one wantw more control over the plot one can plot the individual components. E.g.:

````@example tutorial
plot(det.contacts[1], fillcolor = :red, fillalpha = 1, linecolor = :black, camera = (10, 20), size = (900,900))
plot!(det.contacts[2], fillcolor = :green, fillalpha = 0.25, linecolor = :white)
````

!!! note
    So far, everything is plotted by plotting the individual primitives. Thus, usually, nicer plots are produced
    if the CSG constists only of union of primivites and not differences or intersetions.


## Scalar Potential Plots

## Drift Plots

## Waveform Plots