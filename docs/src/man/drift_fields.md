# Drift Fields

The trajectories of charged particles in vacuum follow the electric field. Following Coulomb's law ($\bm{F} = q \bm{E}$, where $\bm{F}$ corresponds to the force experienced by the particle, $q$ the charge and $\bm{E}$ the electric field), charged particles in vacuum would be continuously accelerated until approaching the speed of light (called ballistic transport). However, inside a material, scattering prevents this constant acceleration and leads to a constant drift velocity 

```math
v_{d} = \mu E,
```
where $v_{d}$ is the drift velocity, $\mu$ the mobility and $E$ the electric field strength.

The scattering with matter reduces the energy of the charge carriers and might also deviate their trajectory from the electric field direction. E.g., in crystals, charges tend to drift along the crystal axes. Since the resulting drift field is not necessarily aligned with the electric field, the electron or hole mobility needs to be expressed as a 3x3 tensor

```math
v_{i} =  \mu_{ij} \cdot E_{j},
```

where $\mu_{ij}$ is the so-called mobility tensor. The mobility varies for different materials and depends also on other parameters such as temperature, impurity concentration and the electric field strenght, as explained later.

Electrons and holes have different mobilities, resulting in different drift fields. There are several models for the mobility tensor of electrons and holes in certain materials. Right now, a drift model for constant charge drift, the electric field charge drift model, and a drift model for high purity germanium crystal, the [ADL/Bruyneel](https://www.sciencedirect.com/science/article/pii/S0168900206015166) model, are implemented in the SSD package. Additionally, the implementation of an own model is possible and simple.

## Electric field charge drift model 

The simplest drift model is the constant charge drift model, which describes the trajectories of electrons and holes moving along the electric field lines. The mobility parameter in this case is a scalar $\pm$ 1 m²/(Vs) ($\+$ for holes, and $\-$ for electrons), and thus, the velocity field has the same (or opposite) direction as the electric field. Even though this model describes a rather ideal and theoretical system, it is useful to obtain and use the electric field vectors as velocity vectors, i.e., to calculate the electric field lines.

To set the electric field charge drift model for the simulation:

```julia
charge_drift_model = ElectricFieldChargeDriftModel()
set_charge_drift_model!(simulation, charge_drift_model) 
calculate_drift_field!(simulation)
```

## ADL charge drift model

The drift velocity of charge carriers inside a germanium crystal can be modelled if the dependency of the mobility parameter on the electric field is known. However, this relation depends as well on the crystal orientation. Germanium detectors have a cubic diamond lattice structure, with <100>, <110> and <111> as principal directions. The dependency of the charge carrier's mobility can be disentagled by parametrizing the mobility along the crystal axes. Assuming that the electric field is aligned to one of the axis, the longitudinal velocity of the charge carrier in the direction $\textit{l}$ can be obtained through the parametrization proposed by [G.F. Knoll](https://www.wiley.com/en-us/Radiation+Detection+and+Measurement%2C+4th+Edition-p-9780470131480):

```math
v_l = \frac{\mu_0 E}{(1 + (E/E_0 )^{\beta})^{1/ \beta})} - \mu_{n} E
```
The parametrization values $\mu_{0}$, $E_{0}$, $\beta$ and $\mu_{n}$ for both electrons and holes have been obtained by [B. Bruyneel](https://www.sciencedirect.com/science/article/pii/S0168900206015166) by measuring the drift velocities in the <100> and <111> directions in a crystal germanium at 78 K. The parameters are written in the config file "drift_velocity_config.json", saved in `<package_directory>/src/ChargeDriftModels/ADL/`.

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

However, a second parametrization is needed to calculate the drift field for configurations where the electric field is not aligned to the principal axes. Different parametrization models have been used for [electrons](https://www.sciencedirect.com/science/article/pii/S0168900299012863) and [holes](https://www.sciencedirect.com/science/article/pii/S0168900206015166), and are implemented in the ADLChargeDriftModel(). In order to perform the calculation of the drift field, a config file with the parametrization values (the "drift_velocity_config.json, with Bruyneel's data, or a user-made one), has to be passed as an argument to the ADLChargeDriftModel function.

```julia
#Commands to find the directory of the config files in the user's computer
ssd_src_path = dirname(Base.find_package("SolidStateDetectors"))
ssd_adl_config_file_path = joinpath(ssd_src_path, "ChargeDriftModels/ADL")
user_adl_charge_drift_model_config_filename = "user_charge_drift_model.json"
path_to_adl_config_file = joinpath(ssd_adl_config_file_path, user_adl_charge_drift_model_config_filename)

#call ADL model with the desired config file 
user_adl_charge_drift_model =ADLChargeDriftModel(path_to_adl_config_file)
calculate_drift_field!((simulation, charge_drift_model) 
```

#####################################################not changed yet######################

While experimental evidence suggests that charge mobility of electrons and holes in the crystal is temperature dependent, the dependency law has not yet been established. Several models have been proposed to reproduce the experimental behaviour, and some examples of them can be found in the same directory `<package_directory>/src/ChargeDriftModels/ADL/`. The examples include a linear model, a Boltzmann model and a power-law model. To use these models in the calculation of the drift field, the name of the config file and the temperature must be given in calling the function. Say we want to use the Boltzmann model at a temperature of 100 K:

```julia
charge_drift_model = ADLChargeDriftModel("drift_velocity_config_boltzmann.json", 100)  
set_charge_drift_model!(simulation, charge_drift_model) 
```


If no temperature is given as a parameter, the calculations will be performed at a default temperature of 80 K.

It should be taken into account that these models have not been proved, and they should only be used as a matter of interpretation.

## User-made charge drift model

User-made charge drift models can also be used in the calculation of the charge velocities. The parameters in the original temperature independent charge drift model (ADL) can be modulated and saved in another config files. The config files must keep the structure and the json format to be interpreted by the function.

The following generic commands allow to obtain the directory where the charge drift models are saved in the user's computer and use them:

```julia
ssd_src_path = dirname(Base.find_package("SolidStateDetectors"))
ssd_adl_config_file_path = joinpath(ssd_src_path, "ChargeDriftModels/ADL")
user_adl_charge_drift_model_config_filename = "user_charge_drift_model.json"
path_to_adl_config_file = joinpath(ssd_adl_config_file_path, user_adl_charge_drift_model_config_filename)

user_adl_charge_drift_model =ADLChargeDriftModel(path_to_adl_config_file)
set_charge_drift_model!(simulation, charge_drift_model) 

```
