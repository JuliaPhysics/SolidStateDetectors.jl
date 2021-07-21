# Weighting Potentials

The weighting potential is a theoretical potential that describes what fraction of a charge at position $\vec{r}$ is seen by a contact, $C_i$. The weighting potential can take values between 0, i.e. the charge is not seen by the electrode, and 1, i.e. the charge is collected by the electrode. 
```math
\nabla \left( \epsilon_r(\vec{r}) \nabla \Phi_i^w(\vec{r})\right) = 0\\
\Phi_i^w(\vec{r})\vert_{C_j} = \left\{ \begin{array}{ll} 1, & \text{if $i = j$} \\ 0, & \text{if $i \neq j$} \end{array} \right.
```
where $\Phi_i^w$ is the electric potential and $\epsilon_r$ is the dielectric distribution.


The net charge induced on each electrode $C_i$, $Q_i$, by electrons and hole with absolute charge $Q$ is given by the Schockley-Ramo theorem:
```math
Q_i = Q \left( \sum\limits_\text{holes} \Phi_i^w(\vec{r}_h) -  \sum\limits_\text{electrons} \Phi_i^w(\vec{r}_e) \right)
```



## Simulation Algorithm

The weighting potential for an electrode is internally calculated with the same function as the electric potential when calling `calculate_weighting_potential!(simulation, contact_id)`. The differences are that $\rho(\vec{r})$ is set to zero and that the boundary conditions of fixed values on the contacts are adapted.