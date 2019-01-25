# Electric Potentials

## Simulation Algorithm

The electric potential is calculated through [successive over relaxation](https://en.wikipedia.org/wiki/Successive_over-relaxation) (SOR).

The calculation is based on Gauss' law in matter
```math
\nabla \epsilon_r(\vec{r}) \nabla \varphi(\vec{r}) = \dfrac{\rho(\vec{r})}{\epsilon_0}\,,
```
which is numerically solved on a 3-dimensional adaptive red/black grid.

The red/black division allows for multithreading and the adaptive grid saves computation time since it only increases the grid point density in areas where it is critical.

For now, only cylindrical coordinates are supported. 

## Multithreading

To use multiple threads for the simulation, the environement variable `JULIA_NUM_THREADS` must be set before Julia is started. In case of bash this is done through
```bash
export JULIA_NUM_THREADS=4
```
Note that the user still has to parse the number of threads to the function as a keyword `nthreads`. See [`calculate_electric_potential`](@ref) and [`calculate_weighting_potential`](@ref).


## Example

To simulate a detector the `detector` must first be defined ([Detectors](@ref)). 

```@example electric_potential
using Plots; pyplot()
using SolidStateDetectors

detector = SolidStateDetector(SSD_examples[:InvertedCoax]);
```

Then, the electric potential can be calculated through the function [`calculate_electric_potential`](@ref). It returns a collection struct [`SolidStateDetectors.PotentialSimulationSetup`](@ref) in which the electric potential is stored. Also it stores the grid, the charge density, the dielectric distribution and the point types [`SolidStateDetectors.PointTypes`](@ref).
The electric potential can be extracted via the function [`SolidStateDetectors.ElectricPotential()`](@ref). The struct ElectricPotential also holds the corresponding grid [`SolidStateDetectors.Grid`](@ref). 

```@example electric_potential
E_pot_setup = calculate_electric_potential(detector);
E_pot = ElectricPotential(E_pot_setup);
```



A plot recipe for the struct [`SolidStateDetectors.ElectricPotential`](@ref) exists so the result can be visualized through
```@example electric_potential
plot(E_pot, θ=0)
```

Since this is a fully θ-symmetric detector, the 3-dimensional grid has only 1 tick in the polar coordinate θ. But due to effects of the crystal axes the grid has to be extended in θ in order to deterime the electric and drift fields. This can be done through [`SolidStateDetectors.extent_2D_grid_to_3D!`](@ref):

```@example electric_potential
E_pot = SolidStateDetectors.get_2π_potential(E_pot, n_points_in_θ = 36);
```

Now it can also be plotted in the r-θ plane: 
```@example electric_potential
plot(E_pot, z=0.04)
```
