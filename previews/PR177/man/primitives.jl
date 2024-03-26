# ## Volume Primitives

using SolidStateDetectors
import SolidStateDetectors.ConstructiveSolidGeometry as CSG
using SolidStateDetectors.ConstructiveSolidGeometry: ClosedPrimitive
using SolidStateDetectors.ConstructiveSolidGeometry: Box
using Plots; gr()
Plots.reset_defaults()
T = Float64;

# ### List of YAML example configuration files for Primitives
# Under `SolidStateDetectors.jl/examples/example_primitive_files` there 
# are some examples how to define the different primitives 
# via the YAML format:

path_to_example_primitives_config_files = joinpath(dirname(dirname(pathof(SolidStateDetectors))), "examples", "example_primitive_files") 
example_primitives_config_filenames = readdir(path_to_example_primitives_config_files) 
for fn in example_primitives_config_filenames
    println(fn)
end

# ### Box
cfn = joinpath(path_to_example_primitives_config_files, "Box.yaml")
print(open(f -> read(f, String), cfn))

# Load the primitive from the configuration file via `CSG.Geometry`
box = CSG.Geometry(T, cfn)
plot(box)

# ### Cone:
# #### Tube
cfn = joinpath(path_to_example_primitives_config_files, "Cone_tube.yaml")
print(open(f -> read(f, String), cfn))

# Load the primitive from the configuration file via `CSG.Geometry`
cone = CSG.Geometry(T, cfn)
plot(cone)

# #### VaryingTube
cfn = joinpath(path_to_example_primitives_config_files, "Cone.yaml")
print(open(f -> read(f, String), cfn))

# Load the primitive from the configuration file via `CSG.Geometry`
cone = CSG.Geometry(T, cfn)
plot(cone)

# ### Ellipsoid
# #### Sphere
cfn = joinpath(path_to_example_primitives_config_files, "Ellipsoid_full_sphere.yaml")
print(open(f -> read(f, String), cfn))

# Load the primitive from the configuration file via `CSG.Geometry`
ellipsoid = CSG.Geometry(T, cfn)
plot(ellipsoid)

# ### Torus
cfn = joinpath(path_to_example_primitives_config_files, "Torus.yaml")
print(open(f -> read(f, String), cfn))

# Load the primitive from the configuration file via `CSG.Geometry`
torus = CSG.Geometry(T, cfn)
plot(torus, zlims = [-6,6], camera = (40, 55))

# ### Prism
# #### Hexagonal Prism
cfn = joinpath(path_to_example_primitives_config_files, "RegularPrism_hexagon.yaml")
print(open(f -> read(f, String), cfn))

# Load the primitive from the configuration file via `CSG.Geometry`
prism = CSG.Geometry(T, cfn)
plot(prism)
