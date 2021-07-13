# Configuration Files

The detector, its surroundings and symmetries can be specified in configuration files.

SSD supports YAML and JSON as formats for the configuration files.

## Example Configuration Files

Several example configuration files can be found under

`<package_directory>/examples/example_config_files/`.

They are accessible through a dictionary, `SSD_examples`, defined in the package:
```@example general
using SolidStateDetectors
keys(SSD_examples) # dictionary holding the full path to the corresponding config JSON files
```

They can be loaded via
```@example general
path_to_config_file = SSD_examples[:InvertedCoax]
sim = Simulation(path_to_config_file)
```