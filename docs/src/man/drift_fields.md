# Drift Fields

Charged particles in vacuum move along the electric field lines under Coulomb's force, $\bm{F} = q \bm{E}$, where $\bm{F}$ corresponds to the force experienced by the particle, $q$ is the charge of the particle and $\bm{E}$ is the electric field. Charged particles in vacuum would be continuously accelerated until approaching the speed of light (called ballistic transport), however, inside a material, scattering prevents this constant acceleration and leads to a constant drift velocity 

```math
v_{d} = \mu E,
```
where $v_{d}$ is the drift velocity, $\mu$ the mobility and $E$ the electric field strength.

The scattering with matter not only limits the absolute drift velocity, it might also deviate the trajectories from the electric field lines: e.g., in crystals, the principal axes orientation has an impact on the resulting drift trajectory. The influence of the scattering on the drift trajectories can be expressed by a 3x3 tensor, the so-called mobility tensor $\mu_{ij}$, which transforms the electric field, $E$, into the drift field, $v_{i}$:

```math
v_{i} =  \mu_{ij} \cdot E_{j}.
```

The mobility varies for different materials and depends also on other parameters such as temperature, impurity concentration and on the electric field strenght, as explained later.

Electrons and holes have different mobilities, resulting in different drift fields. There are several models for the mobility tensor of electrons and holes in certain materials. Right now, two models are implemented. The first one is a pseudo-drift model, the [`ElectricFieldChargeDriftModel`](@ref), which just takes the electric field vectors as drift vectors, see section [Electric Field Charge Drift Model](@ref). The second one, [`ADLChargeDriftModel`](@ref), is a drift model for high purity germanium, see section [ADL Charge Drift Model](@ref). However, the implementation of an own model is possible and explained in section [Custom Charge Drift Model](@ref).

## Electric Field Charge Drift Model

The [`ElectricFieldChargeDriftModel`](@ref) describes a system in which electrons and holes move along the electric field lines. In this case, the mobility is a scalar $\pm$ 1 m²/(Vs) ($+$ for holes, and $-$ for electrons), and thus, the velocity field has the same (or opposite) direction as the electric field. Even though this model does not describe reality, it is useful in some cases to use the electric field vectors as velocity vectors.

In order to set the `ElectricFieldChargeDriftModel` for the simulation, the precision type of the calculation `T` (`Float32` or `Float64`) has to be given as an argument. Note that `T` has to be of the same precision type of the simulation:

```julia
T = SolidStateDetectors.get_precision_type(simulation) # e.g. Float32
charge_drift_model = ElectricFieldChargeDriftModel(T)
set_charge_drift_model!(simulation, charge_drift_model) 
calculate_drift_fields!(simulation)
```

## ADL Charge Drift Model

In high-purity germanium, the mobility cannot be expressed by a simple scalar quantity. Germanium has a cubic diamond lattice structure with $\langle$100$\rangle$, $\langle$110$\rangle$ and $\langle$111$\rangle$ as principal directions. Along these axes, the charge drift is parallel to the electric field. However, the longitudinal drift velocity, $v_{l}$, is not equally fast on the three axes. 

On each axes, $v_{l}$ can be described through the parametrization proposed by  [C. Canali et al.](https://ieeexplore.ieee.org/document/1478102), which was later expanded by [L. Mihailescu](https://www.sciencedirect.com/science/article/pii/S0168900299012863):

```math
v_l = \frac{\mu_0 E}{(1 + (E/E_0 )^{\beta})^{1/ \beta}} - \mu_{n} E.
```

The parameters $\mu_{0}$, $E_{0}$ and $\beta$ differ for electrons and holes and $\mu_{n}$ is only relevant for electrons. These parameters were obtained by [B. Bruyneel](https://www.sciencedirect.com/science/article/pii/S0168900206015166) by measuring the drift velocities of electrons and holes in the $\langle$100$\rangle$ and $\langle$111$\rangle$ directions in high purity germanium at a temperature of 78 K. These parameters are stored in a json config file, "drift\_velocity\_config.json", located in `<package_directory>/src/ChargeDriftModels/ADL/`. The config file is expressed as following:


```json
{
	"phi110": -0.785398,
	"drift": {
		"velocity": {
			"model": "Bruyneel2006",
			"parameters": {
				"e100": {
					"mu0": 3.8609,
					"beta": 0.805,
					"E0": 51100.0,
					"mun": -0.0171
				},
				"e111": {
					"mu0": 3.8536,
					"beta": 0.641,
					"E0": 53800.0,
					"mun": 0.0510
				},
				"h100": {
					"mu0": 6.1824,
					"beta": 0.942,
					"E0": 18500.0
				},
				"h111": {
					"mu0": 6.1215,
					"beta": 0.662,
					"E0": 18200.0
				}
			}
		}
	}
}
```

where the parameters are stored under the keys `e100`, `e111`, `h100` and `h111`, in which `e` and `h` stand for electrons and holes, respectively, and `100` and `111`, for the principal axes $\langle$100$\rangle$ and $\langle$111$\rangle$. Currently, in `SolidStateDetectors.jl` the $\langle$001$\rangle$ axis is fixed to be the Z-axis of the coordinate system of the simulation. The orientation of the crystal is set through the `phi110` parameter, which fixes the angle in radiants between the $\langle$110$\rangle$ principal direction of the crystal and the X-axis.



If the electric field is not aligned with any of the crystal axes, the charge drift velocity is not necessarily aligned with the electric field. In the [`ADLChargeDriftModel`](@ref), two models are implemented to describe the charge drift of electrons and holes between the axes. Detailed information about the charge drift models is provided in the papers from [L. Mihailescu et al. ](https://www.sciencedirect.com/science/article/pii/S0168900299012863) for electrons and from [B.Bruyneel et al.](https://www.sciencedirect.com/science/article/pii/S0168900206015166) for holes.


In order to perform the calculation of the drift field, a json config file containing the parametrization values like the "drift\_velocity\_config.json" (with Bruyneel's data or modified values), has to be passed as an argument to the `ADLChargeDriftModel` function. The precision of the the calculation `T` (`Float32` or `Float64`) has to be given as a keyword `T`. Note that `T` has to be of the same type as the chosen in the simulation:

```julia
T = SolidStateDetectors.get_precision_type(simulation) # e.g. Float32
charge_drift_model = ADLChargeDriftModel("<path_to_ADL_json_config_file>", T=T)
set_charge_drift_model!(simulation, charge_drift_model) 
calculate_drift_fields!(simulation)
```


The values from the default config file correspond to germanium at 78 K. Calculations of the drift field at other temperatures are also supported by the ADL charge drift model. While experimental observations suggest that the charge mobility of electrons and holes in the crystal is temperature dependent, the dependency law has not yet been established. Several models have been proposed to reproduce the experimental behaviour, and some examples of them can be found in the directory `<package_directory>/src/ChargeDriftModels/ADL/`. The examples include a linear model, a Boltzmann model and a power-law model. To use these models in the calculation of the drift fields, the corresponding config file, the temperature and the precision must be given to the function. E.g., in order to use the Boltzmann model at a temperature of 100 K:

```julia
T = SolidStateDetectors.get_precision_type(simulation) # e.g. Float32
charge_drift_model = ADLChargeDriftModel("<path_to_drift_velocity_config_boltzmann.json>", T = T, temperature = 100)  
set_charge_drift_model!(simulation, charge_drift_model) 
calculate_drift_fields!(simulation)
```


If no temperature is given as a parameter, the calculations will be performed at a default temperature of 80 K.

It should be noted that the correct model has not yet been identified, and the parameters inside these config files -besides the default ADL one- are just educated guesses.


## Custom Charge Drift Model

The user can implement and use his own drift model.

The first step is to define a `struct` for the model which is a subtype of `SolidStateDetectors.AbstractChargeDriftModel`:

```julia
using SolidStateDetectors
using SolidStateDetectors: SSDFloat, AbstractChargeDriftModel
using StaticArrays

struct CustomChargeDriftModel{T <: SSDFloat} <: AbstractChargeDriftModel{T} 
    # optional fields to parameterize the model
end
```

The second step is to define two methods (`getVe` for electrons and `getVh` for holes), which perform the transformation of an electric field vector, `fv::SVector{3,T}`, into a velocity vector.
Note, that the vectors are in cartesian coordinates, independent of the coordinate system (cartesian or cylindrical) of the simulation. 

```julia
function SolidStateDetectors.getVe(fv::SVector{3, T}, cdm::CustomChargeDriftModel)::SVector{3, T} where {T <: SSDFloat}
    # arbitrary transformation of fv
    return -fv
end

function SolidStateDetectors.getVh(fv::SVector{3, T}, cdm::CustomChargeDriftModel)::SVector{3, T} where {T <: SSDFloat}
    # arbitrary transformation of fv
    return fv
end
```

Then, one can apply the model to the simulation:
```julia
T = SolidStateDetectors.get_precision_type(simulation) # e.g. Float32
cdm = CustomChargeDriftModel{T}()
set_charge_drift_model!(simulation, cdm)
calculate_drift_fields!(simulation)
```

