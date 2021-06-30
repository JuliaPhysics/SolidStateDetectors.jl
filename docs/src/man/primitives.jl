using SolidStateDetectors
path_to_example_primitives_config_files = joinpath(dirname(dirname(pathof(SolidStateDetectors))), "examples", "example_primitive_files") # hide
example_primitives_config_filenames = readdir(path_to_example_primitives_config_files) # hide
import SolidStateDetectors.ConstructiveSolidGeometry as CSG
using SolidStateDetectors.ConstructiveSolidGeometry: ClosedPrimitive
using SolidStateDetectors.ConstructiveSolidGeometry: Box
using Plots
gr(xlabel = "X", ylabel = "Y", zlabel = "Z"); 
T = Float64;

# # Box
box = Box{T, ClosedPrimitive}(hX = 1.0, hY = 2.0, hZ = 3.0)
plot(box)
#=
YAML configuration file format:
```yaml
box:
  x:
    from: -1.0
    to: 1.0
  y:
    from: -2.0
    to: 2.0
  z:
    from: -3.0
    to: 3.0
```
=#

# # Sphere