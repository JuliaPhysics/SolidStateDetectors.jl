# Capacitances

The mutual capacitance, $c_{ij}$, between two contacts $i$ and $j$ can be calculated via
```julia
calculate_mutual_capacitance(sim, (i, j))
```

The weighting potentials of the two contacts must have already been calculated
as they are needed in the [calculation](https://www.jpier.org/PIERB/pierb92/01.21011501.pdf):
```math
c_{ij} = \epsilon_0 \int_{World} \nabla \Phi_i^w(\vec{r}) ϵ_r(\vec{r}) \nabla \Phi_j^w(\vec{r}) d\vec{r}
```

The function returns the elements of the Maxwell capacitance matrix, $\mathbf{C}$,
which fulfills
```math
\vec{\mathbf{Q}} = \mathbf{Q} \cdot \vec{\mathbf{V}}
```

The hole matrix can also be obtained via 
```julia
calculate_capacitance_matrix(sim)
```