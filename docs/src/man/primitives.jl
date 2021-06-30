using SolidStateDetectors
import SolidStateDetectors.ConstructiveSolidGeometry as CSG
using SolidStateDetectors.ConstructiveSolidGeometry: ClosedPrimitive
using SolidStateDetectors.ConstructiveSolidGeometry: Box
using Plots
T = Float64;

# ## List of YAML example configuration files for Primitives
# Under `SolidStateDetectors.jl/examples/example_primitive_files` there 
# are some examples how to define the different primitives 
# via the YAML format:

path_to_example_primitives_config_files = joinpath(dirname(dirname(pathof(SolidStateDetectors))), "examples", "example_primitive_files") 
example_primitives_config_filenames = readdir(path_to_example_primitives_config_files) 
for fn in example_primitives_config_filenames
    println(fn)
end

# # Box
cfn = joinpath(path_to_example_primitives_config_files, "Box.yaml")
print(open(f -> read(f, String), cfn))

# Load the primitive from the configuration file via `CSG.Geometry`
box = CSG.Geometry(T, cfn)
plot(box)

# # Ellipsoid
cfn = joinpath(path_to_example_primitives_config_files, "Sphere.yaml")
print(open(f -> read(f, String), cfn))

# Load the primitive from the configuration file via `CSG.Geometry`
ellipsoid = CSG.Geometry(T, cfn)
plot(ellipsoid)
