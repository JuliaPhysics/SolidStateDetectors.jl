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