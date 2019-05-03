# Electric Potentials

## Simulation Algorithm

The electric potential is calculated through [successive over relaxation](https://en.wikipedia.org/wiki/Successive_over-relaxation) (SOR).

The calculation is based on Gauss' law in matter

```math
\nabla \epsilon_r(\vec{r}) \nabla \varphi(\vec{r}) = \dfrac{\rho(\vec{r})}{\epsilon_0}\,,
```

which is numerically solved on a 3-dimensional adaptive red/black grid.

The red/black division allows for multithreading and the adaptive grid saves computation time since it only increases the grid point density in areas where it is critical.

Both cylindrical and cartesian coordinates are supported.

## Multithreading

To use multiple threads for the simulation, the environement variable `JULIA_NUM_THREADS` must be set before Julia is started. In case of bash this is done through

```bash
export JULIA_NUM_THREADS=4
```

Note that the user still has to parse the number of threads to the function as a keyword `nthreads`. See [`calculate_electric_potential!`](@ref) and [`calculate_weighting_potential!`](@ref).

## Example

First the Simulation Object must be initialized from a config file. This is the central structure that will be gradually extended with all fields and results. Initially it only holds the holds the information about the geometry and electric properties of all objects in the 'world'. This information is summarized in the field 'detector'

```@example electric_potential
using Plots; pyplot()
using SolidStateDetectors

mySSD = Simulation(SSD_examples[:InvertedCoax]);
detector = mySSD.detector
```


The electric potential can be calculated through the function [`calculate_electric_potential!`](@ref). It fills the respective field of the `Simulation`. Also it stores the grid, the charge density, the dielectric distribution and the point types [`PointTypes`](@ref).
The electric potential can be accessed via [`mySSD.electric_potential`](@ref). The electric potential is a struct itself and also holds the corresponding grid [`Grid`](@ref).

```@example electric_potential
calculate_electric_potential!(mySSD);
E_pot = mySSD.electric_potential;
```

A plot recipe for the struct [`ElectricPotential`](@ref) exists so the result can be visualized through

```@example electric_potential
plot(E_pot, φ=0)
```

It can also be plotted in the r-φ plane:

```@example electric_potential
plot(E_pot, z=0.04)
```
