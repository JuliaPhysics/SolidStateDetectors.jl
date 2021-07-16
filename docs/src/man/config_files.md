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