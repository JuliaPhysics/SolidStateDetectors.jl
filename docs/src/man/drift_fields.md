# Drift Fields

The drift trajectories of charged particles in vacuum follow the electric field. Following Coulomb's law ($\bm{F} = q \bm{E}$, where $\bm{F}$ corresponds to the force experienced by the particle, $q$ the charge and $\bm{E}$ the electric field), charged particles in vacuum would be continuously accelerated until approaching the speed of light (called ballistic transport). However, inside a material, scattering prevents this constant acceleration and leads to a constant drift velocity 

```math
v_{d} = \mu E,
```
where $v_{d}$ is the drift velocity, $\mu$ the mobility and $E$ the electric field strength.

The scattering with matter reduces the energy of the charge carriers and in some cases might also deviate their trajectory from the electric field direction: e.g., in crystals, the principal axes orientation influence the resulting drit velocity. Since the resulting drift field is not necessarily aligned with the electric field, the electron or hole mobility needs to be expressed as a 3x3 tensor

```math
v_{i} =  \mu_{ij} \cdot E_{j},
```

where $\mu_{ij}$ is the so-called mobility tensor. The mobility varies for different materials and depends also on other parameters such as temperature, impurity concentration and on the electric field strenght, as explained later.

Electrons and holes have different mobilities, resulting in different drift fields. There are several models for the mobility tensor of electrons and holes in certain materials. Right now, a drift model for constant charge drift -the so-called [electric field charge drift model](../api.md#SolidStateDetectors.ElectricFieldChargeDriftModel)-, and a drift model for high purity germanium crystal, the [ADL charge drift model](../api.md#SolidStateDetectors.ADLChargeDriftModel), are implemented. Additionally, the implementation of an own model is possible and simple.

## Electric field charge drift model 

The simplest drift model is the constant charge drift model, which describes a system in which electrons and holes move exclusively along the electric field lines. The mobility parameter in this case is a scalar $\pm$ 1 m²/(Vs) ($+$ for holes, and $-$ for electrons), and thus, the velocity field has the same (or opposite) direction as the electric field. Even though this model describes a rather ideal and theoretical system, it is useful in some cases to obtain and use the electric field vectors as velocity vectors.

In order to set the electric field charge drift model for the simulation, the precision of the calculation `T` (Float32 or Float64) has to be given as an argument. Note that `T` has to be of the same type as in the simulation:

```julia
T = Float32
charge_drift_model = ElectricFieldChargeDriftModel(T)
set_charge_drift_model!(simulation, charge_drift_model) 
calculate_drift_fields!(simulation)
```

## ADL charge drift model

The drift velocity of charge carriers inside a germanium crystal can be modelled if the dependencies of the mobility parameter on external parameters are known. Particularly in the case of the dependency on the electric field, this relation depends on its turn on the crystal orientation. Germanium detectors have a cubic diamond lattice structure, with $\langle$100$\rangle$, $\langle$110$\rangle$ and $\langle$111$\rangle$ as principal directions. The dependency of the charge carriers mobility can be disentangled by parametrizing the mobility along the crystal axes. Assuming that the electric field is aligned with one of the axes, the longitudinal velocity of the charge carrier in the direction $\textit{l}$ can be obtained through the parametrization proposed by [C. Canali et al.](https://ieeexplore.ieee.org/document/1478102) and later expanded by [L. Mihailescu](https://www.sciencedirect.com/science/article/pii/S0168900299012863):

```math
v_l = \frac{\mu_0 E}{(1 + (E/E_0 )^{\beta})^{1/ \beta})} - \mu_{n} E.
```
The parametrization values $\mu_{0}$, $E_{0}$ and $\beta$ for both electrons and holes, and $\mu_{n}$ only for electrons, have been obtained by [B. Bruyneel](https://www.sciencedirect.com/science/article/pii/S0168900206015166) by measuring the drift velocities in the $\langle$100$\rangle$ and $\langle$111$\rangle$ directions in germanium at 78 K. These parameters are written in the config file "drift\_velocity\_config.json", located in `<package_directory>/src/ChargeDriftModels/ADL/`. The config file is expressed as following:

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

OPTIONAL = (where the parameters are stored under the keys `e100`, `e111`, `h100` and `h111`, in which the _e_ and _h_ stand for electrons and holes, respectively, and 100 and 111, for the principal axes $\langle$100$\rangle$ and $\langle$111$\rangle$.) While in the ADL charge drift model the $\langle$001$\rangle$ principal direction of the crystal is fixed on the Z-axis, the crystal can rotate in the x-y plane. The orientation of the crystal is set through the `phi110` parameter, which fixes the angle in radiants between the $\langle$110$\rangle$ principal direction of the crystal and the X-axis.


A second parametrization is needed to calculate the drift field for points in which the electric field is not aligned to the principal axes. Different parametrization models have been used for [electrons](https://www.sciencedirect.com/science/article/pii/S0168900299012863) and [holes](https://www.sciencedirect.com/science/article/pii/S0168900206015166), and are implemented in the [`ADLChargeDriftModel`](../api.md#SolidStateDetectors.ADLChargeDriftModel). In order to perform the calculation of the drift field, a json config file with the parametrization values (the "drift\_velocity\_config.json, with Bruyneel's data, or a custom one), has to be passed as an argument to the ADLChargeDriftModel function. The precision of the the calculation `T` (Float32 or Float64) has to be given as a keyword `T`. Note that `T` has to be of the same type as the chosen in the simulation:

```julia
T = Float32
charge_drift_model = ADLChargeDriftModel("<path_to_ADL_json_config_file>", T=T)
set_charge_drift_model!(simulation, charge_drift_model) 
calculate_drift_fields!(simulation)
```


The values from the default config file correspond to germanium at 78 K. Calculations of the drift field at other temperatures are also supported by the ADL charge drift model. While experimental observations suggest that charge mobility of electrons and holes in the crystal is temperature dependent, the dependency law has not yet been established. Several models have been proposed to reproduce the experimental behaviour, and some examples of them can be found in the same directory `<package_directory>/src/ChargeDriftModels/ADL/`. The examples include a linear model, a Boltzmann model and a power-law model. To use these models in the calculation of the drift fields, the corresponding config file, the temperature and the precision must be given in calling the function. Say we want to use the Boltzmann model at a temperature of 100 K:

```julia
T = Float32
charge_drift_model = ADLChargeDriftModel("<path_to_drift_velocity_config_boltzmann.json>", T = T, temperature = 100)  
set_charge_drift_model!(simulation, charge_drift_model) 
calculate_drift_fields!(simulation)
```


If no temperature is given as a parameter, the calculations will be performed at a default temperature of 80 K.

It should be noted that the correct model has not yet been identified, and the parameters inside these config files -besides the default ADL one- are just educated guesses.


#####################################################not changed yet######################

## Custom charge drift model

The user caa implement and use his own drift model.

to be done...
