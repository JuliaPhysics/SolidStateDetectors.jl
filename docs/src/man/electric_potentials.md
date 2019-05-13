# Electric Potentials

## Simulation Algorithm

The electric potential is calculated through [successive over relaxation](https://en.wikipedia.org/wiki/Successive_over-relaxation) (SOR).

The calculation is based on Gauss' law in matter

```math
\nabla \left( \epsilon_r(\vec{r}) \nabla \Phi(\vec{r})\right) = \dfrac{\rho(\vec{r})}{\epsilon_0}\,,
```
where $\Phi$ is the electric potential, $\rho$ is the charge density,
$\epsilon_r$ is the dielectric distribution and $\epsilon_0$ is the dielectric constant of the vacuum.

The equation is numerically solved on a 3-dimensional adaptive red/black grid.

The red/black division allows for multithreading and the adaptive grid saves computation time since it only increases the grid point density in areas where it is critical.

Cylindrical and cartesian coordinates are supported.

## Multithreading

To use multiple threads for the simulation, the environement variable `JULIA_NUM_THREADS` must be set before Julia is started. In case of bash this is done through

```bash
export JULIA_NUM_THREADS=4
```

Note that the user still has to parse the number of threads to the function as a keyword `nthreads`.

ToDo...