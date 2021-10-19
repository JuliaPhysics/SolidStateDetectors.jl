# Plotting

In order to plot geometries or simulation results, the user will need to load the [Plots.jl](https://github.com/JuliaPlots/Plots.jl) package.

````@example tutorial
using Plots
using SolidStateDetectors
using Unitful

T = Float32
sim = Simulation{T}(SSD_examples[:InvertedCoax]);
simulate!(sim, convergence_limit = 1e-6, refinement_limits = [0.2, 0.1, 0.05, 0.01]);
````

Besides the usual [Plots.jl](https://github.com/JuliaPlots/Plots.jl) keywords settings, like `size = (900, 900)`,
there are some additional settings to tune the plots which are described in the following sections.

## Detector Plots

The geometry of a detector, together with its environment, can be simply plotted via
````@example tutorial
det = sim.detector
plot(det, size = (900,900))
````

### Optional keywords:
* `seriestype`: Can be `:mesh3d` (default).
* `linewidth`: Sets the line width of the edges of the mesh.
* `linecolor`: Sets the line color of all edges of the mesh.
* `fillcolor`: Sets the face color of all faces.
* `fillalpha`: Sets the alpha value of all faces of the mesh.
* `show_semiconductor`: Whether also the primitives of the semiconductor should be plotted. Default is `false`.
* `show_passives`: Whether also the primitives of the surrounding objects of the detector should be plotted. Default is `true`.

#### How does the plot recipe work?

The detector consists out of a semiconductor, `det.semiconductor`, its contacts, `det.contacts`,
and its surrounding objects, `det.passives`.
By default, see optional keywords, only the contacts and the surrounding objects are plotted
and the color of the contact primitives is defined internally through their `id`.

If one wants more control over the plot one can plot the individual components, e.g.
````@example tutorial
plot(det.contacts[1], fillcolor = :red, fillalpha = 1, linecolor = :black, camera = (10, 20), size = (900,900))
plot!(det.contacts[2], fillcolor = :green, fillalpha = 0.25, linecolor = :white)
````

!!! note
    So far, everything is plotted by plotting the individual primitives. Thus, usually, nicer plots are produced
    if the geometry consists only of union of primitives and not differences or intersections.


## Scalar Potential Plots

When calculating the electric potential, a set of scalar potentials are stored in the `Simulation` object, e.g. `sim.electric_potential`, `sim.point_types`, and `sim.q_eff_imp`. SolidStateDetectors.jl has plot recipes for these scalar potentials.

````@example tutorial
plot(sim.electric_potential)
````

There are several keyword arguments that can be passed to adjust the potential plots.
First, the simulated potentials are simulated on a three-dimensional `Grid`, but the plots will always be a two-dimensional cross section (for a `CylindricalGrid` at constant `r`, `φ` or `z`, for a `CartesianGrid3D` at constant `x`, `y` or `z`).
The cross section can be specified in the `plot` command. A cross section of the `ElectricPotential` at `z = 20mm` can be plotted via
````@example tutorial
plot(sim.electric_potential, z = 20u"mm")
````

In addition to all plot attributes that are implemented in [Plots.jl](https://github.com/JuliaPlots/Plots.jl), plots of scalar potentials can take two additional keyword arguments, `contours_equal_potential` and `full_det`, which both are of type `Bool`:
* `contours_equal_potential`: If this set to `true`, the plot will additionally display contour lines at points with equal potential value. 
* `full_det`: By default, cross sections in `φ` are only displayed for positive radii `r`. If `full_det` is set to `true`, cross sections in `φ` are also extended to "negative" `r`, by additionally plotting the potential values at `φ + 180°` on the left side of the plot.


````@example tutorial
plot(sim.electric_potential, φ = 30u"°", contours_equal_potential = true, full_det = true)
````


## Electric Field Plots

The magnitude of the electric field can be plotted using the same syntax as [Scalar Potential Plots](@ref).

````@example tutorial
plot(sim.electric_field, full_det = true)
````

In addition, SolidStateDetectors.jl offers the possibilities to plot electric field lines via `plot_electric_fieldlines`.
This is done by spawning charges close to the surface of the contacts and simulating their drift
parallel to the electric field, using [`ElectricFieldChargeDriftModel`](@ref).

````@example tutorial
plot_electric_fieldlines!(sim, full_det = true)
````

In addition to all plot attributes that are implemented in [Plots.jl](https://github.com/JuliaPlots/Plots.jl),
the syntax for specifying the cross section and `full_det` (see [Scalar Potential Plots](@ref)),
the plot recipe for electric field lines can additionally understand the following arguments to tune the plot.
* `sampling`: Specifies the steps at which the contacts are sampled to generate equally spaced charges at the surface. The default is `2u"mm"`, but the optimal value depends on the geometry of the detector and contacts. If no unit is given, `sampling` is interpreted in units of meter.
* `offset`: The charges are created on the surface and have to be moved slightly inside the semiconductor to be able to drift (charges that are in the contacts will not drift). This keyword defines how much the charges will be moved inside along the normal vector of the surface. The default is `0.5u"mm"`, but the optimal value again depends on the detector geometry. Unitless quantities are interpreted in units of meter.
* `spacing`: Fine sampling can cause the plot two be too cluttered. The `spacing` keyword allows to skip some field lines. The default is `1`, which means that every line is plotted. If only every second line is to be plotted, set `spacing = 2`, etc.
* `skip_contact`: Detectors will usually have positively and negatively biased contacts. The charges should only be spawned on either of them, but not both. This keyword defines which contact should be skipped in the charge spawn algorithm. Default is `1`, i.e. that no charges will be spawned at the surface of contact with `id = 1`.
* `max_nsteps`: After all, the generation of electric field lines is based on the charge drift code, which requires a maximum number of steps, which is set to `5000` by default here. If the drift ends before the charges reach a contact, the electric field line will not be fully displayed. In this case, it is recommended to increase `max_nsteps`.



The electric field can be plotted using the . The electric field strength is plotted using `plot(sim.electric_field)`, whereas the electric field lines can be plotted on top of that using  `plot_electric_fieldlines!(sim)`.

Minimum working example:
````@example tutorial
plot(sim.electric_field, full_det = true)
plot_electric_fieldlines!(sim, full_det = true)
````

## Event Plots

SolidStateDetectors.jl also provides plot recipes to display charge drift paths resulting from simultaneous energy deposits
in the semiconductor body of the detector. To simulate and plot an individual event, the `Event` struct is recommended.

````@example tutorial
evt = Event([CartesianPoint{T}(0.01,0.01,0.075)])
simulate!(evt, sim)
plot(sim.detector, size = (500,500), label = "")
plot!(evt.drift_paths)
````

The drift path plots can be modified using the keywords implemented in [Plots.jl](https://github.com/JuliaPlots/Plots.jl):
````@example tutorial
plot(sim.detector, size = (500,500), label = "")
plot!(evt.drift_paths, linewidth = 2, linestyle = :dash, markersize = 6)
````

There are also plot recipes for plotting the simulated waveforms of the `Event`:
````@example tutorial
plot(evt.waveforms, unitformat = :slash, label = "Contact ".*string.((1:2)'), legend = :topleft)
````

The length of the waveforms is given by the length of the charge drift.
By default, no baseline and no tail are added to the waveforms. However, this might be desired in waveforms plots.
The waveforms can be extended by calling [`add_baseline_and_extend_tail`](@ref) on the waveforms:
````@example tutorial
plot(add_baseline_and_extend_tail.(evt.waveforms,0,400), linewidth = 8, linestyle = :dash, label = "", unitformat = :slash)
````
