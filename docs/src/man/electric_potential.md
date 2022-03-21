# Electric Potential

The electric potential is given by Gauss' law in matter
```math
\nabla \left( \epsilon_r(\vec{r}) \nabla \Phi(\vec{r})\right) = - \dfrac{\rho(\vec{r})}{\epsilon_0}\,,
```
where $\Phi$ is the electric potential, $\rho$ is the charge density,
$\epsilon_r$ is the dielectric distribution and $\epsilon_0$ is the dielectric constant of the vacuum.


## Simulation Algorithm

The electric potential is calculated through [successive over relaxation](https://en.wikipedia.org/wiki/Successive_over-relaxation) when calling `calculate_electric_potential!(sim)`. The equation is numerically solved on a three-dimensional adaptive red/black grid. The red/black division allows for multithreading and the adaptive grid saves computation time since it only increases the grid point density in areas where it is critical. To use multiple threads for the simulation, the environment variable `JULIA_NUM_THREADS` must be set before Julia is started. In case of bash this is done through
```bash
export JULIA_NUM_THREADS=4
```

!!! note "GPU"
    The electric potential can be calculated on GPUs. See [GPU Support in Field Calculations](@ref).


At the beginning of the simulation, each grid point, $(i,j,k)$ is assigned its dielectric constant, $\epsilon_r(\vec{r}_{i,j,k})$, as well as its effective charge, $Q_\text{eff} = \rho(\vec{r}_{i,j,k}) \cdot V_{i,j,k} / \epsilon_0$, where $V_{i,j,k}$ is the volume assigned to the grid point $(i,j,k)$.

These quantities are stored in the fields `sim.q_eff_imp` and `sim.ϵ_r` and can be plotted using
```julia
using SolidStateDetectors
using Plots 
sim = Simulation(SSD_examples[:InvertedCoax])
apply_initial_state!(sim, ElectricPotential)
plot(
  plot(sim.q_eff_imp),
  plot(sim.ϵ_r)
)
```


## Impurity Densities

One contribution to the charge density $\rho(\vec{r})$ is the impurity density of the semiconductor of a detector.
Some simple impurity density profiles are already implemented in SolidStateDetectors.jl and can be easily accessed in the configuration files. Note that all impurity densities are given in units of **atoms / particles** per volume.

!!! note "Sign of the impurity density"
    The sign of the impurity density determines whether the semiconductor is p-type or n-type.

    p-type region <-> negative sign: Holes are the majority carriers and are free to move and diffuse into the n-type region. 
        Electrons are fixed in the lattice. Thus, a negative fixed space charge density is left behind in the depleted p-type region.

    n-type region <-> positive sign: Electrons are the majority carriers and are free to move and diffuse into the p-type region. 
        Holes are fixed in the lattice. Thus, a positive fixed space charge density is left behind in the depleted n-type region.


### Constant Impurity Density

A constant impurity density throughout the detector volume can be modeled with `ConstantImpurityDensity`.
In the configuration files, `constant` impurity densities are defined with the `value` of the constant impurity density, i.e.
```yaml
impurity_density:
  name: constant
  value: 1e10cm^-3 # => 10¹⁹ m⁻³
```
If no units are given, `value` is interpreted in units of `units.length`$^{-3}$.
They are converted to SI units (m$^{-3}$) internally.


### Linear Impurity Density

An impurity density with a linear gradient can be modeled with `LinearImpurityDensity`.
In the configuration files, `linear` impurity densities are defined with an `init` (initial) value and `gradient` along
each Cartesian direction (`x`, `y` and `z`), e.g.
```yaml
impurity_density:
  name: linear
  z:
    init: 1e10cm^-3
    gradient: 1e10cm^-4
```
or
```yaml
impurity_density:
  name: linear
  x:
    init: 0
    gradient: 1e10
  z:
    init: 0
    gradient: 1e10
```
In the first example, the `init` value corresponds to the value at `z = 0` whereas the gradient points towards positive `z`.
In the second example, the impurity density is 0 at the origin of the coordinate system, whereas the gradient of the impurity density profile points in $\langle101\rangle$ direction.
If no units are given, `init` is parsed in units of `units.length`$^{-3}$ and `gradient` in units of `units.length`$^{-4}$.


### Cylindrical Impurity Density

An impurity density with a radial gradient can be modeled with `CylindricalImpurityDensity`.
In the configuration files, `cylindrical` impurity densities are defined with an `init` (initial) value and `gradient` along
each cylindrical spatial direction (`r` and `z`), e.g.
```yaml 
impurity_density:
  name: cylindrical
  r:
    init: 1e10cm^-3
    gradient: 1e10cm^-4
```
Here, the impurity density at the origin is $10^{10}$cm$^{-3}$ and it increases radially with the gradient $10^{10}$cm$^{-4}$.
If no units are given, `init` is parsed in units of `units.length`$^{-3}$ and `gradient` in units of `units.length`$^{-4}$.


### Custom Impurity Density

The source code for the previously introduced impurity densities can be found [here](https://github.com/JuliaPhysics/SolidStateDetectors.jl/blob/master/src/ImpurityDensities). More complex impurity density profiles can be defined by the user.
Each custom impurity density is a new `struct` and subtype of `SolidStateDetectors.AbstractImpurityDensity`
and needs a method `SolidStateDetectors.get_impurity_density` that returns the impurity density at a given point `pt`.


#### Example 1: Radially Oscillating Impurity Density

```julia
using SolidStateDetectors: AbstractChargeDensity, CartesianVector, AbstractCoordinatePoint
import SolidStateDetectors: get_impurity_density

# new struct for translated impurity densities
struct OscillatingImpurityDensity{T} <: AbstractImpurityDensity{T}
    wavelength::T 
    amplitude::T
    offset::T
end

# add get_charge_density for the newly defined charge density model
function SolidStateDetectors.get_impurity_density(ocdm::OscillatingImpurityDensity{T}, pt::AbstractCoordinatePoint{T})::T where {T}
    cyl_pt = CylindricalPoint(pt) # convert point to a CylindricalPoint
    return ocdm.offset + ocdm.amplitude * sin(2π * cyl_pt.r / ocdm.wavelength)
end
```

#### Example 2: Translating Existing Impurity Densities

```julia
using SolidStateDetectors: AbstractImpurityDensity, CartesianVector, AbstractCoordinatePoint
import SolidStateDetectors: get_impurity_density

# new struct for translated impurity densities
struct TranslatedImpurityDensity{T} <: AbstractImpurityDensity{T}
    impurity_density_model::AbstractImpurityDensity{T}
    translate::CartesianVector{T}
end

# add get_impurity_density for the newly defined impurity density model
function SolidStateDetectors.get_impurity_density(tcdm::TranslatedImpurityDensity{T}, pt::AbstractCoordinatePoint{T})::T where {T}
    translated_pt::CartesianPoint{T} = CartesianPoint(pt) - tcdm.translate
    return get_impurity_density(tcdm.impurity_density_model, translated_pt)
end
```

#### Example 3: P-N Junction

Have a look at [Advanced Example: Custom Impurity Profile](@ref).

## Charge Densities

Another contribution to the charge density $\rho(\vec{r})$ can be charged surfaces or volumes that can be modeled using passive objects. The same profiles as for impurity densities are defined here that can be accessed similarly, i.e.
```yaml
charge_density:
  name: constant
  value: 1e-10 # => 10⁻¹⁰ C/m⁻³
```
Note that, in contrast to impurity densities, charge densities are given in units of the **elementary charge** per volume.