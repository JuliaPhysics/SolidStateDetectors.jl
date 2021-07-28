# Electric Field

The electric field is calculated from the electric potential for each grid point ``(i,j,k)`` when calling `calculate_electric_field!(sim)`. The field vector components on each grid point are the means of the electric field in each direction calculated as finite differences.

On a Cartesian grid, this results in:
```math
\vec{E}^{i,j,k} = \left( \mathcal{E}_x^{i,j,k}, \mathcal{E}_y^{i,j,k}, \mathcal{E}_z^{i,j,k} \right)^{\mathsf{T}}~\hspace{10pt},
```
where
```math
\begin{aligned}
	\mathcal{E}_x^{i,j,k} &= \dfrac{1}{2}\left( \dfrac{\Phi_{i+1,j,k}-\Phi_{i+1,j,k}}{x_{i+1} - x_{i}} + \dfrac{\Phi_{i,j,k}-\Phi_{i-1,j,k}}{x_{i} - x_{i-1}} \right)\hspace{10pt},\\
	\mathcal{E}_y^{i,j,k} &= \dfrac{1}{2}\left( \dfrac{\Phi_{i,j+1,k}-\Phi_{i,j,k}}{y_{j+1} - y_{j}} + \dfrac{\Phi_{i,j,k}-\Phi_{i,j-1,k}}{y_{j} - y_{j-1}} \right)\hspace{10pt},\\
	\mathcal{E}_z^{i,j,k} &= \dfrac{1}{2}\left( \dfrac{\Phi_{i,j,k+1}-\Phi_{i,j,k}}{z_{k+1} - z_{k}} + \dfrac{\Phi_{i,j,k}-\Phi_{i,j,k-1}}{z_{k} - z_{k-1}} \right)\hspace{10pt}.
\end{aligned}
```


On a cylindrical grid, the calculation is:
```math
\vec{E}^{i,j,k} = \left( \mathcal{E}_r^{i,j,k}, \mathcal{E}_{\varphi}^{i,j,k}, \mathcal{E}_z^{i,j,k} \right)^{\mathsf{T}}~\hspace{10pt},
```
where
```math
\begin{aligned}
	\mathcal{E}_r^{i,j,k} &= \dfrac{1}{2}\left(\dfrac{\Phi_{i+1,j,k}-\Phi_{i,j,k}}{r_{i+1} - r_{i}} + \dfrac{\Phi_{i,j,k}-\Phi_{i-1,j,k}}{r_{i} - r_{i-1}}\right)\hspace{10pt},\\
	\mathcal{E}_{\varphi}^{i,j,k} &= \dfrac{1}{2 r_{i}}\left(\dfrac{\Phi_{i,j+1,k}-\Phi_{i,j,k}}{\varphi_{j+1} - \varphi_{j}} + \dfrac{\Phi_{i,j,k}-\Phi_{i,j-1,k}}{\varphi_{j} - \varphi_{j-1}}\right)\hspace{10pt},\\
	\mathcal{E}_z^{i,j,k} &= \dfrac{1}{2}\left( \dfrac{\Phi_{i,j,k+1}-\Phi_{i,j,k}}{z_{k+1} - z_{k}} + \dfrac{\Phi_{i,j,k}-\Phi_{i,j,k-1}}{z_{k} - z_{k-1}} \right)\hspace{10pt}.
\end{aligned}
```

This discrete electric field is interpolated (via [Interpolations.jl](https://github.com/JuliaMath/Interpolations.jl)) 
during the drift in order to get the electric field at the current position of the charge carrier.

The calculated electric field is stored as a field in the Simulation object, i.e. `sim.electric_field`, consisting of the actual field vectors for each grid point stored in `sim.electric_field.data` and the corresponding grid in `sim.electric_field.grid`. 

## Visualization

The electric field can be plotted using the [Plots.jl](https://github.com/JuliaPlots/Plots.jl) package. The electric field strength is plotted using `plot(sim.electric_field)`, whereas the electric field lines can be plotted on top of that using  `plot_electric_fieldlines!(sim)`.

Minimum working example:
```julia
using SolidStateDetectors
using Plots
sim = Simulation(SSD_examples[:InvertedCoax])
calculate_electric_potential!(sim)
calculate_electric_field!(sim)
plot(sim.electric_field, full_det = true)
plot_electric_fieldlines!(sim, full_det = true)
```