# Plotting

In order to plot geometries or simulation results, the user will need to load the [Plots.jl](https://github.com/JuliaPlots/Plots.jl) package.

First, we have to calculate something in order to demonstrate the plotting tools:
````@setup tutorial
using Plots
using SolidStateDetectors
using Unitful
gr(fmt = :png) #hide

T = Float32
sim = Simulation{T}(SSD_examples[:InvertedCoax]);
simulate!(sim, convergence_limit = 1e-6, refinement_limits = [0.2, 0.1, 0.05, 0.01]);
````
````julia
using Plots
using SolidStateDetectors
using Unitful

T = Float32
sim = Simulation{T}(SSD_examples[:InvertedCoax]);
simulate!(sim, convergence_limit = 1e-6, refinement_limits = [0.2, 0.1, 0.05, 0.01]);
````

Besides the usual [Plots.jl](https://github.com/JuliaPlots/Plots.jl) keywords settings, like `size = (500, 500)`,
there are some additional settings to tune the plots which are described in the following sections.

## Detector Plots

The geometry of a detector, together with its environment, can be simply plotted via
````@example tutorial
det = sim.detector
plot(det, size = (500, 500))
````
### Plot Styles

The style of the detector plot can be controlled via the `seriestype` keyword. There are 4 available styles:

````@example tutorial
plot(
      plot(det, seriestype = :csg, title = ":csg"),
      plot(det, seriestype = :wireframe, title = ":wireframe"),
      plot(det, seriestype = :mesh3d, title = ":mesh3d"), 
      plot(
            det, seriestype = :samplesurface, n_samples = 100, 
            markersize = 2, markeralpha = 0.03, title = ":samplesurface"
            ), 
      layout = (1,4), size = (800,200), legend = false, ticks = false, 
      guide = "", zlims = (-0.005,0.1)
)
````

The `seriestype` can be set to `:csg` (default), `:wireframe`, `:mesh3d`, or `:samplesurface`. `:csg` plots a wireframe on top of a mesh (with no mesh gridlines). For fastest plotting use either `:wireframe` or `:mesh3d` and consider changing `n_arc` (see [Optional Keywords](@ref)). `:csg`, `:wireframe`, and `:mesh3d` are all mesh-based. For geometries containing differences or intersections the recommended seriestype is `:samplesurface`. This can be seen by plotting the detector's semiconductor below. `:samplesurface` is marker-based. The marker density is set low by default for speed. For increased plot fidelity use `n_samples = 100` and `markersize = 2`.

````@example tutorial
plot(
      plot(det.semiconductor, seriestype = :csg, title = ":csg"),
      plot(det.semiconductor, seriestype = :wireframe, title = ":wireframe"),
      plot(det.semiconductor, seriestype = :mesh3d, title = ":mesh3d"), 
      plot(
            det.semiconductor, seriestype = :samplesurface, n_samples = 100, 
            markersize = 2, markeralpha = 0.03, title = ":samplesurface"
            ), 
      layout = (1,4), size = (800,200), legend = false, ticks = false, 
      guide = "", zlims = (-0.005,0.1)
)
````

!!! note
    So far, when using mesh-based seriestypes (`:csg`, `:wireframe`, `:mesh3d`), plots are produced by plotting whole primitives. Thus, usually, nicer plots are produced if the geometry consists only of unions of primitives and not differences or intersections. If the geometry contains differences, the resulting negative geometries are plotted with thinner wireframe lines and/or with semi-transparent white mesh faces depending on the `seriestype` used. 
    
### How does the plot recipe work?

The detector consists of a semiconductor, `det.semiconductor`, its contacts, `det.contacts`,
and its surrounding objects, `det.passives`.
By default (see [Optional Keywords](@ref)), only the contacts and the surrounding objects are plotted
and the color of the contact primitives is defined internally through their `id`.

Components can also be plotted individually for enhanced style handling. Additionally, the units of the axes are set by calling a `plot` command with units beforehand.
````@example tutorial
plot(u"cm", u"cm", u"cm")
plot!(det.semiconductor, st = :samplesurface, n_samples = 100, markersize = 2,
      camera = (40, 55), size = (500, 500))
plot!(det.contacts[1], st = :mesh3d, linewidth = 0.5, fillcolor = :white)
plot!(det.contacts[2], st = :wireframe, n_vert_lines = 5)
````

### Optional Keywords
* `show_semiconductor`: Will display the semiconductor if set to `true`. Default is `false`.
* `show_passives`: Will display the objects surrounding the detector if set to `true`. Default is `true`.
* `seriestype`: Can be `:csg` (default), `:wireframe`, `:mesh3d`, or `:samplesurface`.
* `linewidth`: Sets the line width of the edges of the mesh gridlines when using `seriestype = :mesh3d`. When using `seriestype = :csg` or `seriestype = :wireframe`, `linewidth` sets the line width of the wireframe.
* `linecolor`: Sets the line color of the edges of the mesh gridlines when using `seriestype = :mesh3d`. When using `seriestype = :csg` or `seriestype = :wireframe`, `linecolor` sets the line color of the wireframe.
* `fillcolor`: Sets the face color of all faces of the mesh.
* `fillalpha`: Sets the alpha value of all faces of the mesh.
* `markercolor`: Sets the marker color.
* `markersize`: Sets the marker size. `markersize = 4` is the default value.
* `markeralpha`: Sets the alpha value for markers.
* `n_arc`: Controls the discretization of curved objects in a mesh. Each full ellipse is divided into `n_arc` segments. Partial ellipses are drawn with a proportional number (`n_arc*f` with `f<1`) of segments.`n_arc = 40` is the default value. Smaller `n_arc` values will result in faster plotting, specially if the geometry contains tori or ellipsoids.
* `n_vert_lines`: Controls the number of wireframe "vertical" lines in a mesh. `n_vert_lines = 2` is the default value. A maximum of `n_arc*f` vertical lines can be drawn on each curved object. This keyword is ignored by polygons.
* `n_samples`: Controls the marker density. `n_samples = 40` is the default value. Reduce this value for faster plotting. Consider increasing `markersize` and/or `markeralpha` to compensate for the visual impact of lower marker densities. Note that the marker density is intended to be even across all dimensions. Therefore, visual distortions will occur if the aspect ratio of the axes is far from unity.

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
plot( sim.electric_potential, φ = 30u"°", 
      contours_equal_potential = true, full_det = true, 
      linecolor = :white, levels = 34)
````


## Electric Field Plots

The magnitude of the electric field can be plotted using the same syntax as [Scalar Potential Plots](@ref).

````@example tutorial
plot(sim.electric_field, full_det = true, clims = (0, 5e5))
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
* `skip_contact`: Detectors will usually have positively and negatively biased contacts. The charges should only be spawned on either of them, but not both. This keyword defines which contact should be skipped in the charge spawn algorithm. Default is `1`, i.e. that no charges will be spawned at the surface of contact with `id = 1`.
* `max_nsteps`: After all, the generation of electric field lines is based on the charge drift code, which requires a maximum number of steps, which is set to `5000` by default here. If the drift ends before the charges reach a contact, the electric field line will not be fully displayed. In this case, it is recommended to increase `max_nsteps`.



In case the electric field line plots do not look good, adjust `sampling`, `offset` and `max_nsteps` until obtaining the desired result.
````@example tutorial
plot(sim.electric_field, full_det = true, size = (500,500), clims = (0,5e5))
plot_electric_fieldlines!(sim, full_det = true, sampling = 3u"mm", offset = 2u"mm")
````

## Event Plots

SolidStateDetectors.jl also provides plot recipes to display charge drift paths resulting from simultaneous energy deposits
in the semiconductor body of the detector. To simulate and plot an individual event, the `Event` struct is recommended.

````@example tutorial
evt = Event([CartesianPoint{T}(0.01,0.01,0.075)], [2u"MeV"])
simulate!(evt, sim)
plot(sim.detector, size = (500,500), label = "")
plot!(evt.drift_paths)
````

The drift path plots can be modified using the keywords implemented in [Plots.jl](https://github.com/JuliaPlots/Plots.jl).
The units of the axes can be set by calling a `plot` command with units beforehand:
````@example tutorial
plot(u"mm", u"mm", u"mm")
plot!(sim.detector, size = (500,500), label = "")
plot!(evt.drift_paths, linewidth = 2, linestyle = :dash, markersize = 6)
````

There are also plot recipes for plotting the simulated waveforms of the `Event`:
````@example tutorial
plot(evt.waveforms, unitformat = :slash, label = "Contact ".*string.((1:2)'), legend = :topleft)
````

The length of the waveforms is given by the length of the charge drift.
By default, no baseline and no tail are added to the waveforms. However, this might be desired in waveforms plots.
The waveforms can be extended by calling [`add_baseline_and_extend_tail`](@ref) on the waveforms.
Again, the units of the axes can be set by calling a `plot` command with units before plotting the waveforms.
````@example tutorial
plot(u"µs", u"fC")
plot!(add_baseline_and_extend_tail.(evt.waveforms,0,400), 
      linewidth = 4, linestyle = :dash, 
      label = "", unitformat = :slash)
````
