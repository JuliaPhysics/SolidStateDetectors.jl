# Grids

The electric potential is calculated via [successive over relaxation](https://en.wikipedia.org/wiki/Successive_over-relaxation)
on a 3-dimensional grid. SolidStateDetectors.jl can calculate the electric potential on a Cartesian or a cylindrical grid:
* `CartesianGrid`: 3 axes: `x`-, `y`- and `z`-axis.
* `CylindricalGrid`: 3 axes: `r`-, `φ`- and `z`-axis.

The system to simulate, e.g. a cryostat with a detector inside, is called "the world".
This world is divided into a set of discrete points, called the `grid::Grid`.
It is defined through three axes: `grid.axes`. Each of the three axes is divided into a
discrete number of points (ticks): $N_1, N_2, N_3$.  
The linear combinations of those points form the set of all
$N_\mathrm{gp} = N_1 \times N_2 \times N_3$ grid points.

Each axis is defined as a [`SolidStateDetectors.DiscreteAxis`](@ref).

## DiscreteAxis

A [`SolidStateDetectors.DiscreteAxis`](@ref) defines the axis of a dimension of a grid, e.g. the `x`-axis of a `CartesianGrid`. 
A `DiscreteAxis`, e.g. `ax::DiscreteAxis`, stores the boundary conditions (reflecting / periodic)
at the endpoints of the axis as well as the discrete ticks, `ax.ticks`, of the axis in the interval `ax.interval`.
The interval also specifies whether the endpoints of the interval are included in the ticks or not
via `:closed` or `:open` in the type of the interval.
The ticks do not need to be evenly spaced, allowing for an adaptive refinement of the grid
in areas where the gradient of the electric potential is large.


## Grid Boundary Conditions

In the potential calculation, the potential value at a given grid point is determined by the potential values of its nearest neighbors.
At the boundaries of the world, there are no nearest neighbors.
This requires choosing boundary conditions on how to continue the potential beyond the volume of the world.

````@example grid
using Plots, LaTeXStrings #hide
xs = [1,2,3,4,5,6,7,8,9] #hide
ys = [0.38,0.46,0.48,0.44,0.35,0.28,0.2,0.15,0.1] #hide
nothing #hide
````

In order to explain the different boundary conditions, assume a `DiscreteAxis` with `N` ticks, `x1` to `xN`, where `x1` and `xN` define the boundary of the world in that dimension. In this example, `N = 9`:
````@example grid
vline([extrema(xs)...], color = :black, ls = :dash, label = "") #hide
vline!([extrema(xs) .+ (-1,1)...], color = :darkgray, ls = :dash, label = "") #hide
scatter!(xs, ys, label = "Potential", color = 1) #hide
P = plot!(size = (600,200), xlims = (-0.5,10.5), ylims = (0,0.55), xlabel = "Grid ticks", ylabel = "Potential values", xticks = (1:9, latexstring.("x\$_".*string.(1:9).*"\$")), xtickfontsize = 12, margin = 5Plots.mm, ygrid = false, legend = (0.73,1)) #hide
````

For the potential calculation, extended ticks are added left and right of the outermost grid ticks. 
In this example, they are labeled as `x0` and `xN+1`.
The choice of boundary condition influences which potential values these virtual extended ticks acquire in the potential calculation.

`periodic`: This assumes the potential to be periodic beyond the world volume, i.e. `ϕ(x0) = ϕ(xN)` and `ϕ(xN+1) = ϕ(x1)`:
````@example grid
A = deepcopy(P) #hide
plot!(A, xticks = (0:10, [L"\"x$_0$\""; latexstring.("x\$_".*string.(1:9).*"\$"); L"\"x$_{10}$\""])) #hide
plot!(A, [extrema(xs) .+ (0,1)...], [ys[begin],ys[begin]], arrow = arrow(:both, :closed), label = "", color = :orange) #hide
plot!(A, [reverse(extrema(xs)) .- (0,1)...], [ys[end],ys[end]], arrow = arrow(:both, :closed), label = "", color = :orange) #hide
scatter!(A, [0,length(xs)+1], [ys[end], ys[begin]], label = "periodic", color = :orange) #hide
````

`reflecting`: This assumes mirror-symmetry in the potential with respect to the world boundary, i.e. `ϕ(x0) = ϕ(x2)` and `ϕ(xN+1) = ϕ(xN-1)`:
````@example grid
A = deepcopy(P) #hide
plot!(A, xticks = (0:10, [L"\"x$_0$\""; latexstring.("x\$_".*string.(1:9).*"\$"); L"\"x$_{10}$\""])) #hide
plot!(A, [minimum(xs) .+ (-1,1)...], [ys[begin+1],ys[begin+1]], arrow = arrow(:both, :closed), label = "", color = :purple) #hide
plot!(A, [maximum(xs) .+ (-1,1)...], [ys[end-1],ys[end-1]], arrow = arrow(:both, :closed), label = "", color = :purple) #hide
scatter!([0,length(xs)+1], [ys[begin+1], ys[end-1]], label = "reflecting", color = :purple) #hide
````

`fixed` (or `fix`): This sets the potential beyond the world volume to the potential at the world boundary, i.e. `ϕ(x0) = ϕ(x1)` and `ϕ(xN+1) = ϕ(xN)`:
````@example grid
A = deepcopy(P) #hide
plot!(A, xticks = (0:10, [L"\"x$_0$\""; latexstring.("x\$_".*string.(1:9).*"\$"); L"\"x$_{10}$\""])) #hide
plot!(A, [minimum(xs) .+ (-1,0)...], [ys[begin],ys[begin]], arrow = arrow(:both, :closed), label = "", color = :green) #hide
plot!(A, [maximum(xs) .+ (1,0)...], [ys[end],ys[end]], arrow = arrow(:both, :closed), label = "", color = :green) #hide
scatter!(A, [0,length(xs)+1], [ys[begin], ys[end]], label = "fix / fixed", color = :green) #hide
````

`infinite` (or `inf`): Assumes the potential to follow a relationship beyond the world volume, where the ratio between neighboring ticks is expected to stay constant, i.e. `ϕ(x0)/ϕ(x1) = ϕ(x1)/ϕ(x2)` and `ϕ(xN+1)/ϕ(xN) = ϕ(xN)/ϕ(xN-1)`:
````@example grid
A = deepcopy(P) #hide
plot!(A, xticks = (0:10, [L"\"x$_0$\""; latexstring.("x\$_".*string.(1:9).*"\$"); L"\"x$_{10}$\""])) #hide
plot!(A, -1:0.01:11, x -> ys[end] * (ys[end]/ys[end-1])^(x - xs[end]), color = :red, ls = :dash, label = "") #hide
plot!(A, -1:0.01:11, x -> ys[begin] * (ys[begin+1]/ys[begin])^(x - xs[begin]), color = :red, ls = :dash, label = "") #hide
scatter!(A, [0,length(xs)+1], [ys[begin]^2/ys[begin+1], ys[end]^2 / ys[end-1]], label = "inf / infinite", color = :red) #hide
````



## Grid Initialization

The grid can be specified in the configuration files, see [Grid](@ref), of `sim::Simulation`.

There is a constructor method for grid: [`Grid(sim::Simulation)`](@ref).

In [`calculate_electric_potential!`](@ref),[`calculate_weighting_potential!`](@ref) and [`simulate!(sim::Simulation)`](@ref)
the calculation of the potentials starts with an initial grid, which can be passed to these functions via the keyword `grid`.
If no grid is passed, the grid is generated via [`Grid(::Simulation)`](@ref).

The two keywords `max_tick_distance` and `max_distance_ratio` can also be passed to the functions
[`calculate_electric_potential!`](@ref),[`calculate_weighting_potential!`](@ref) and [`simulate!(sim::Simulation)`](@ref)
where they are internally forwarded.

### Initialization of the Grid

1) Important points of the objects of the simulation are obtained, e.g. the corners of a `Box`-primitive.
    The coordinates of these points are used to generate the ticks of each axis of the grid.
2) For each axis, additional ticks are added based on the keyword `max_distance_ratio::Real = 5`:
    If the ratio between the distances of an axis tick to its neighboring ticks is larger than `max_distance_ratio`
    (or smaller than `inv(max_distance_ratio)`), additional ticks are inserted such that in the end
    all those ratios are within the interval [`inv(max_distance_ratio)`, `max_distance_ratio`].
3) Finally, additional ticks are added for each axis if the distance between two ticks is larger
    than a threshold distance specified via the keyword `max_tick_distance`.

`max_tick_distance` can either be a `Quantity`, e.g. `1u"mm"`, or a Tuple of Quantities, e.g. `(1u"mm", 0.2u"cm", 3u"mm")`,
to set it for each axis of the grid separately. If `max_tick_distance` is `missing`, one fourth of the axis length is used.

See also [`Grid(::Simulation)`](@ref).

## Grid Refinement

After the potential has converged to equilibrium on the initial grid, the grid can be refined (multiple times).

The refinement can be tuned via the keyword `refinement_limits` in
[`calculate_electric_potential!`](@ref),[`calculate_weighting_potential!`](@ref) and [`simulate!(sim::Simulation)`](@ref).

It defines the maximum (relative to applied bias voltage) allowed differences
of the potential values of neighbored grid points in each dimension for each refinement.
It can be specified in different ways, see e.g. [`calculate_electric_potential!`](@ref).

One simple example would be `refinement_limits = [0.2, 0.1, 0.05]`.
This would mean that the grid would be refined three times and the refinement limit would be
the same for each dimension of the grid in each refinement.
In the first refinement, the refinement limit would be 0.2. Thus, if the bias voltage of the detector in the simulation
is `1000 V`, the maximum allowed potential difference between two grid points would be `200 V`.
Let's say there is a potential difference of `500 V` between the two grid points at `(i,j,k)` and `(i,j+1,k)`.
Then, two (`floor(Int, 500 / 200) = 2`) ticks are added in the second axis between the
previous axis ticks, `grid.axes[2].ticks[j]` and `grid.axes[2].ticks[j+1]`,
which results in $2 \times N_1 \times N_3$ new additional grid points.
The potential values at the added grid points are determined through linear interpolation.
Then, the potential values of all grid points are updated (through the SOR) until convergence
is reached again and the next refinement with `0.1` is executed.

Another keyword can be used to set a minimum allowed distance between two ticks: `min_tick_distance`, see e.g. [`calculate_electric_potential!`](@ref), which prohibits the insertion of new ticks if the new resulting distances between the ticks would be below this limit.

!!! note
    It is usually favourable to make more refinements with smaller differences in the limits.
    Thus, for example, `refinement_limits = [0.2, 0.1, 0.05, 0.03, 0.01]` is usually better than
    `refinement_limits = [0.2, 0.01]`.
