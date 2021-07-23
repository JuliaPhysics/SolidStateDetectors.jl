# Grids

The electric potential is calculated via [successive over relaxation](https://en.wikipedia.org/wiki/Successive_over-relaxation)
on a 3-dimensional grid.

The grid, e.g. `grid::Grid`, does not need to be evenly spaced.
It is defined through the axes of the grid: `grid.axes`,
which are defined as [`SolidStateDetectors.DiscreteAxis`](@ref).

The number of total grid points is the product of the number of ticks of each axis.

SSD can calculate the electric potential on a `CartesianGrid` or a `CylindricalGrid`.

* `CartesianGrid`: 3 axes: x-, y- and z-axis.
* `CylindricalGrid`: 3 axes: r-, φ- and z-axis.

## DiscreteAxis

A [`SolidStateDetectors.DiscreteAxis`](@ref) defines the axis of a dimension of a grid.
E.g. the x-axis of a CartesianGrid.
A `DiscreteAxis`, e.g. `ax::DiscreteAxis`, stores the boundaries conditions (reflecting / periodic)
at the ends of the axis as well as the discrete ticks, `ax.ticks`, of the axis in the interval `ax.interval`.
The interval also specifies whether the endpoints of the interval are included in the ticks or not
via `:closed` or `:open` in the type of the interval.

## Grid Initialization

The grid can be specified in the configuration files, see [Grid](@ref), of `simulation::Simulation`.

There is a constructor method for grid: [`Grid(sim::Simulation)`](@ref).

In [`calculate_electric_potential!`](@ref),[`calculate_weighting_potential!`](@ref) and [`simulate!(sim::Simulation)`](@ref)
the calculation of the fields start with an initial grid, which can be passed to these functions via the keyword `grid`.
If no grid is passed, the grid is generated via [`Grid(::Simulation)`](@ref)`.

The two keywords `max_tick_distance` and `max_distance_ratio` can also be passed to the functions
[`calculate_electric_potential!`](@ref),[`calculate_weighting_potential!`](@ref) and [`simulate!(sim::Simulation)`](@ref)
where they are internally forwarded.
### How the Grid is generated:

1) Important points of the objects of the simulation are obtained. E.g. the corners of a `Box`-primitive.
    The coordinates of these points are used to generate the ticks of each axis of the grid.
2) For each axis additional ticks are added based on the keyword `max_distance_ratio::Real = 5`:
    If the ratio between the distances of an axis tick to its neighboring ticks is large than `max_distance_ratio`
    (or smaller than `inv(max_distance_ratio)`) additional ticks are inserted such that in the end
    all those ratios are within the interval [`inv(max_distance_ratio)`, `max_distance_ratio`].
3) Finally, for each axis additional ticks are added if the distance between to ticks is larger
    than a threshold distance specified via the keyword `max_tick_distance`.

Look at [`Grid(::Simulation)`](@ref) for more details.
## Grid Refinement

After the field has converged on the initial grid, the grid can be refined (multiple times).

The refinement can be tuned via the keyword `refinement_limits` in
[`calculate_electric_potential!`](@ref),[`calculate_weighting_potential!`](@ref) and [`simulate!(sim::Simulation)`](@ref).

It defines the maximum relative (to applied bias voltage) allowed differences
of the potential value of neighbored grid points in each dimension for each refinement.
It can be specified in different ways, see e.g. [`calculate_electric_potential!`](@ref).

One simple example would be `refinement_levels = [0.2, 0.1, 0.05]`.
This would mean that the grid would be refined three times and the refinement limit would be the same for each dimension of the grid.
In the first refinement, the refinement limit would be 0.2. Thus, if the bias voltage of the detector in the simulation
is `1000 V`, the maximum allowed potential difference between to grid points would be `200 V`.
Let's say there is an potential difference of `500 V` between to grid points `(i,j,k)` and `(i,j+1,k)`.
Than, two (`floor(500 V / 200 V)`) ticks would be added in the second axis between the previous axis ticks, `ax.ticks[j]` and `ax.ticks[j+1]`.
The potential values at the added grid points would be determined through linear interpolation.
Than, the potential is updated (trough the SOR) until convergence is reached again.

Another keyword can be used to set a minimum allowed distance between to ticks: `min_tick_distance`, see e.g. [`calculate_electric_potential!`](@ref),
which prohibits the insertion of new ticks if the new resulting distances between the ticks would be below this limit.