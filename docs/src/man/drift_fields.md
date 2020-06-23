# Drift Fields

The trajectories of charged particles in vacuum follow the electric field lines. The electric field is calculated according to Poisson's law through the gradient of the electic potential ($\bm{E} =  -\nabla \phi$). Following Coulomb's law ($\bm{F} = q \bm{E}$), charged particles in vacuum would be accelerated to an ever-increasing velocity (called ballistic transport). However, inside a material, the continual scatter between the electrons or holes with other charged particles (or other structures such as phonons), prevents a constant acceleration and leads to a constant drift velocity for the charged particle

```math
v_{d} = \mu E,
```
where $v_{d}$ is the drift velocity ($m/s$ in IS), $\mu$ the mobility ($m^{2}$/ (V $ \cdot$ s )) and $E$ the magnitude of the electric field (V/m).

The scatter deviates the trajectory of the electrons and holes from the electric field direction. In the case of crystals, their geometry plays a major rol, and the charges tend to drift along the crystal axes. Therefore a parametrization of the electric field on the crystal axes is needed. This parametrization is defined in the charge mobility. Hence, the mobility parameter that relates the velocity field (or drift field) and the electric field line is no longer a scalar parameter, but a 3x3 tensor:

```math
v_{i} =  \mu_{ij} \cdot E_{j},
```

where $\mu$ is the so-called mobility tensor. 

Electrons and holes have different mobilities, and so are their drift fields. There are several models for the mobility tensor of electrons and holes in a crystal; the one implemented in the current package is the [ADL charge drift model](https://www.sciencedirect.com/science/article/pii/S0168900206015166).

## Vacuum charge drift model 
As a matter of an example, a vacuum charge drift model can be selected for the calculation of the electron trajectories in vacuum. The mobility parameter in vacuum is a scalar $(\pm)$ 1 $m^{2}$/ (V $ \cdot$ s )) ($(+) for holes, and $(-) for electrons), thus, the velocity field has the same (or opposite) direction as the electric field.

To set the vacuum charge drift model for the simulation:

```julia
charge_drift_model = VacuumChargeDriftModel()
set_charge_drift_model!(simulation, charge_drift_model) 
apply_charge_drift_model!(simulation)
```
## ADL charge drift model

Inside a crystal, more elaborate charge drift model is required in order to simulate the trajectories of charges.  As mentioned above, the charge drift model implemented in the SSD package is the [ADL charge drift model](https://www.sciencedirect.com/science/article/pii/S0168900206015166). By declaring it, the temperature independent charge drift model is used in the calculations. This model expresses the parameters in a crystal at 78 K.
 

```julia
charge_drift_model = ADLChargeDriftModel() #same as charge_drift_model = ADLChargeDriftModel("drift_velocity_config.json")
set_charge_drift_model!(simulation, charge_drift_model) 
```
The model is written in a json config file, in the directory `<package_directory>/src/ChargeDriftModels/ADL/drift_velocity_config.json`. The parameters are expressed as following:

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
				
} 	"E0": 51100.0,
					"mun": -0.0171
				},
				"e111": {
					"mu0": 3.8536,
					"beta": 0.641,
					"E0": 53800.0,
					"mun": 0.0510
				},
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
		},
		"chrystaltemp": 100
	}
```

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
