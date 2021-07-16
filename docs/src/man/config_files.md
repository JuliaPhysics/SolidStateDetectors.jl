# Configuration Files

The detector, its surroundings and symmetries can be specified in configuration files.

SolidStateDetectors.jl supports YAML and JSON as formats for the configuration files.

## Example Configuration Files

Several example configuration files can be found under

`<package_directory>/examples/example_config_files/`.

They are accessible through a dictionary, `SSD_examples`, defined in the package:
```@example general
using SolidStateDetectors
keys(SSD_examples) # dictionary holding the full path to the corresponding configuration files
```

They can be loaded via
```@example general
path_to_config_file = SSD_examples[:InvertedCoax]
sim = Simulation(path_to_config_file)
```

## General Structure

The configuration files need a minimum of information in order to define the detector, its surroundings and symmetries.

This is a minimum working example of a simple true coaxial detector with two contacts:
```yaml
name: Simple True Coax # optional
units:
  length: mm
  angle: deg
grid:
  coordinates: cylindrical
  axes:
    r: 45
    z:
      from: -40
      to: 40
medium: vacuum
detectors:
- bulk:
    material: HPGe
    geometry:
      tube:
        r: 
          from: 0.5cm
          to: 4cm
        h: 6cm
  contacts:
  - material: HPGe
    channel: 1
    potential: 6000
    geometry:
      tube:
        r:
          from: 5
          to: 5
        h: 60
  - material: HPGe
    channel: 2
    potential: 0
    geometry:
      tube:
        r:
          from: 40
          to: 40
        h: 60
```
It will be used to guide through the different parts of the configuration file.



### Units

Internally, SolidStateDetectors.jl performs its calculations in SI units. However, configuration files can be written in custom units.

The field `units` denotes the standard units with which values will be parsed. Standard units can be defined for `length`, `angle`, `potential` and `temperature`.

In the example above, 
```yaml
units:
  length: mm
  angle: deg
```
will lead to all `length` values to be parsed in units of `mm`, while all `angle` values will be parsed in units of `deg` (degree).

The configuration files also allow for directly passing units to the values that will be parsed using `uparse` from the Unitful.jl package, e.g.
```yaml
units: 
  length: mm
  # ....
tube:
  r: 
    from: 5
    to: 40
  h: 60
```
is equivalent to
```yaml
tube:
  r: 
    from: 5mm
    to: 40mm
  h: 60mm
```
or
```yaml
tube:
  r: 
    from: 0.5cm
    to: 4cm
  h: 6cm
```
In the last example, even if the `length` unit was set to `mm`, the values will be parsed in units of `cm`. Please note to not leave a white space between the value and the unit and to use the Unitful.jl notation.



### Grid 

The calculations are performed on a finite world. To define the world, SolidStateDetectors.jl requires the properties of the grid, which are the coordinate system type and the dimensions. These are defined in the `grid` section of the configuration file, e.g.
```yaml
grid:
 coordinates: cartesian
 axes:
   x: 
     from: -40
     to: 40
   y:
     from: -40
     to: 40
   z:
     from: -40
     to: 40
```

The `coordinates` of the `grid` can be:
- `cartesian` (with the axes `x`, `y` and `z`)
- `cylindrical` (with the axes `r`, `phi` and `z`).

The `axes` field is used to define the dimensions of each axis and, optionally, the boundary handling.
In the example above, the `x`, `y` and `z` axes range from `-40` to `40` units.


#### Grid boundary handling

Symmetries of the world can be used to reduce the calculation only to a fraction of the world. These can be passed as `boundaries` to the different `axes`.

For linear axes (`x`, `y`, `z`), the `boundaries` can be chosen `infinite`, `periodic`, `reflecting`, or `fixed`.

For radial axes (`r`), the boundaries can be chosen `r0`.
If no boundaries are given, the default is `r0` for the left boundary and `infinite` for the right boundary.

For angular axes (`phi`), the boundaries can be chosen `reflecting` or `periodic`. If no boundaries are given, the default is `periodic` for both edges.

All $\varphi$-symmetric configurations can be calculated in 2D if `phi` ranges from `0` to `0` with `periodic` boundary handling, i.e.
```yaml
grid:
 coordinates: cylindrical
 axes:
   r: #...
   phi:
     from: 0
     to: 0
     boundaries: periodic
   z: #...
``` 

All $\varphi$-periodic configurations can be calculated on the fraction of the full $2\pi$ interval, i.e. for a `120°`-periodic system
```yaml
grid:
 coordinates: cylindrical
 axes:
   r: #...
   phi:
     from: 0°
     to: 120°
     boundaries: periodic
   z: #...
``` 

Different boundary handling can be chosen for the `left` and `right` end of the interval, i.e. for a `60°`-periodic system with mirror symmetry at `30°`
```yaml
grid:
 coordinates: cylindrical
 axes:
   r: #...
   phi:
     from: 0°
     to: 30°
     boundaries:
       left: periodic
       right: reflecting
   z: #...
``` 



### Detector Constituents

The detectors for the simulation are defined in an array `detectors`, where each entry corresponds to one detector.
Each detector consists of exactly one `bulk`, a minimum of two `contacts` and, optionally, `passives` and `virtual_drift_volumes`.

```yaml
detectors:
  - name: "Detector 1"
    bulk: #...
    contacts: 
      - # Contact 1
      - # Contact 2 
    passives: 
      - # Passive 1 (optional)
    virtual_drift_volumes:
      - # Virtual Drift Volume 1 (optional)
  - name: "Detector 2"
    bulk: #...
    contacts: 
      - # Contact 1
      - # Contact 2
```

#### Bulk 

An example definition of the bulk looks like this:
```yaml 
bulk:
  material: HPGe
  temperature: 78
  impurity_density: # ...
  charge_drift_model: # ...
  geometry: # ...
```

The different fields of the bulk are:
- `material`: the material of the bulk. This is important to know the electric properties of the bulk for the electric potential calculation. Possible choices are `HPGe` (high-purity germanium) and `Si` (silicon).
- `temperature` (optional): the temperature of the bulk. If no `temperature` is given, the default is 78K for germanium and 293K for all other materials.
- `impurity_density` (optional): the distribution of impurities in the semiconductor material. This has a strong impact on the electric potential calculation. If no `impurity_density` is given, the default is an impurity-free material ($\rho(\vec{r}) = 0$).
- `charge_drift_model` (optional): a model to describe the drift of charge carriers in the semiconductor material. If no `charge_drift_model` is given, the default is `ElectricFieldChargeDriftModel`. Find a detailed description on how to define a charge drift model [here](Custom-Charge-Drift-Model).
- `geometry`: the geometry of the semiconductor object. Find a detailed description on how to define geometries in the section CSG.


## Splitting Configuration Files

Configuration files for complex geometries can get quite long. SolidStateDetectors.jl allows for splitting configuration
files into smaller ones and loading them using the `include` keyword. This feature supports YAML and JSON files.

When including a separate file, the user has to add its file path in the main configuration file at the place it supposed to be added. To identify the file, set the key of this entry to `include`. Here, the user can also give an array of file paths. The file paths can be relative to the path of the configuration file or absolute. When including nested files and using relative paths, please always refer to the last parent file.

Including one file:
```yaml
include : "file_to_be_included.yaml"
```

Including a list of files:
```yaml
include: 
  - "first_file_to_be_included.yaml"
  - "second_file_to_be_included.yaml"
```

Add files to an array in the main configuration file
```yaml
detectors:
  - include: "first_file_in_array.yaml"
  - include: "second_file_in_array.yaml"
  - include: "thrid_file_in_array.yaml"
```

A fully working example can be seen in `SSD_examples[:InvertedCoaxInCryostat]`. Here, the channels, the geometry and other parts are split into separate configuration files.