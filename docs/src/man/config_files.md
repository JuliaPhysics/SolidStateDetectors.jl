# Config Files

## Example Detector Config Files

Currently, there are four predefined example detectors. 

Three of them are cylindrical detectors:

    * Coaxial detector (Coax)
    * Inverted coax detector (IVC)
    * BEGe type detector (BEGe)
The fourth one is in cartesian coordinates:

    * simple cube detector (CGD)

They are all specified in their JSON config files, which can be found under:
`<package_directory>/examples/example_detector_config_files/<config_filename>.json`.

In Julia, their path is already saved in the SolidStateDetectors.jl package in a dictionary:
```julia
    using SolidStateDetectors
    SSD_examples # dictionary holding the full path to the corresponding config JSON files
    SSD_exmpless[:Coax]
```
The keys are: `:Coax`, `:InvertedCoax`, `:BEGe`, `:CGD`.

### Example 1) Inverted Coax

Example minimum config file for an Inverted Coax detector (IVC) plus explanations.
**Remember**, comments are not allowed in JSON files and have to be deleted if you want to use it.
```json
{
    name = "Example Detector",
    ToDo...
}
```

## UserConfigs: Define your own geometry

ToDo...