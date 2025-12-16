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
  value: 1e10cm^-3 # => 10¹⁶ m⁻³
```
If no units are given, `value` is interpreted in units of `units.length`$^{-3}$.
They are converted to SI units (m$^{-3}$) internally.


### Linear Impurity Density

An impurity density with a linear gradient can be modeled with `LinearImpurityDensity`.
In the configuration files, `linear` impurity densities are defined with an `offset` value and `gradient` along
each Cartesian direction (`x`, `y` and `z`), e.g.
```yaml
impurity_density:
  name: linear
  offset: 1e10cm^-3
  gradient: 
    z: 1e10cm^-4
```
or
```yaml
impurity_density:
  name: linear
  offset: 0 # optional
  gradient:
    x: 1e10
    z: 1e10
```
In the first example, the `offset` value corresponds to the impurity density value at the origin of the coordinate system, whereas the gradient points towards positive `z`.
In the second example, the impurity density is 0 at the origin of the coordinate system, whereas the gradient of the impurity density profile points in $\langle101\rangle$ direction.
If no units are given, `offset` is parsed in units of `units.length`$^{-3}$ and `gradient` in units of `units.length`$^{-4}$.


### Cylindrical Impurity Density

An impurity density with a radial gradient can be modeled with `CylindricalImpurityDensity`.
In the configuration files, `cylindrical` impurity densities are defined with an `offset` value and `gradient` along
each cylindrical spatial direction (`r` and `z`), e.g.
```yaml 
impurity_density:
  name: cylindrical
  offset: 1e10cm^-3
  gradient:
    r: 1e10cm^-4
```
Here, the impurity density at the origin is $10^{10}$cm$^{-3}$ and it increases radially with the gradient $10^{10}$cm$^{-4}$.
If no units are given, `offset` is parsed in units of `units.length`$^{-3}$ and `gradient` in units of `units.length`$^{-4}$.

### P-type PN junction Impurity Density

A PN junction impurity model based on lithium thermal diffusion and custom bulk impurity density. The surface lithium density is at the saturation level.
ref: [Dai _et al._ (2023)](https://doi.org/10.1016/j.apradiso.2022.110638)

This density model can be defined in the config file, e.g.
```yaml
impurity_density:
  name: PtypePNjunction
  lithium_annealing_temperature: 623K
  lithium_annealing_time: 18minute
  doped_contact_id: 2
  bulk_impurity_density:
    name: constant
    value: -1e10cm^-3
```
The `impurity_density` needs:
- `name`: the name of the impurity density model, which in this case is `PtypePNjunction`.
- `lithium_annealing_temperature` (optional): lithium annealing temperature, when the lithium is diffused into the crystal. The default value is 623 K.
- `lithium_annealing_time` (optional): lithium annealing time. The default value is 18 minutes.
- `doped_contact_id`:  the doped contact id.
- `bulk_impurity_density`: the density profile of the p-type impurities.

The example above defines a detector with a constant p-type impurity density of
$10^{10}$cm$^{-3}$ and a highly-doped Lithium contact close to the contact with ID 2, created using a lithium annealing temperature of 623 K and an annealing time of 18 minutes.

It is based on the internal implementation of `distance_to_surface`, which might not work for all combinations of primitives.

Users can also define their own method for calculating the distance to the contact in code, e.g. using
```julia
using SolidStateDetectors
using SolidStateDetectors: ConstantImpurityDensity
T = Float64
sim = Simulation{T}(SSD_examples[:TrueCoaxial])

distance_to_contact = pt -> 0.01 - hypot(pt[1], pt[2])
lithium_annealing_temperature::T = 623 # Kelvin
lithium_annealing_time::T = 18*60 # seconds
p_type_density::T = -1e16 # m^-3
bulk_imp_model = ConstantImpurityDensity{T}(p_type_density)

pn_junction_impurity_density = PtypePNJunctionImpurityDensity{T}(lithium_annealing_temperature, lithium_annealing_time, nothing, bulk_imp_model, distance_to_contact)
sim.detector = SolidStateDetector(sim.detector, pn_junction_impurity_density)
```

### Boule Impurity Densities

Boule impurity densities are meant to be used when the impurity density is defined in the boule coordinates, where the z-axis is aligned with the boule growth direction. 
Different models are provided. In each the field `det_z0` is the z-coordinate of the detector origin in boule coordinates. 
The z-direction of the detector is opposite to the z-direction of the boule coordinates.
In this matter the detector impurities are automatically determined from those of the boule, depending on `det_z0`.
```yaml 
impurity_density:
  name: linear_exponential_boule
  a: -1e10cm^-3
  b: -1e9cm^-4
  n: -2e9cm^-3
  l: 5cm
  m: 3cm
  det_z0: 120cm
```
In this example, the impurity density is modeled along the boule as a linear plus an exponential term: $a + b*z + n*e^{(z - l)/m}$ with $z$ in boule coordinates. 
The impurity density of the detector will in turn be modeled by: $a + b*(z_0 - z) + n*e^{(z_0 - z - l)/m}$ with $z$ in detector coordinates.
If no units are given, `a` and `n` are parsed in units of `units.length`$^{-3}$, `b` in units of `units.length`$^{-4}$ and `l`, `m` and `det_z0` in units of `units.length`.

### Correcting Impurity Densities

When simulating real detectors, the simulated depletion voltage often differs from the measured value. This discrepancy arises primarily from the high uncertainty in impurity density measurements. A straightforward approach to improve the agreement is to apply a scaling factor and an offset to the impurity density. 
The scaling factor accounts for potential systematic errors in impurity density measurements. 
The offset is motivated by the thermal release of impurities from deep hole traps, which leads to a shift in the effective impurity density throughout the detector after biasing. In steady-state operation, this effect can be approximated as a uniform offset
These terms can be added to the configuration file under `corrections`. 
```yaml 
impurity_density:
  name: linear_exponential_boule
  a: -1e10cm^-3
  b: -1e9cm^-4
  n: -2e9cm^-3
  l: 5cm
  m: 3cm
  det_z0: 120cm
  corrections:
    scale: 0.9
    offset: -1e9cm^-3
```
If no units are given, `offset` is parsed in units of `units.length`$^{-3}$.
Corrections can be applied to any impurity density model. Given an originally calculated impurity density $\rho$, the corrected impurity density used in the simulation is $f*\rho+t$ where $f$ and $t$ are the scale and offset respectively.

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
  value: 1e-10 # => 10⁻¹⁰ C/m³
```
Note that, in contrast to impurity densities, charge densities are given in units of the **Coulomb** per volume (instead of elementary charge per volume).