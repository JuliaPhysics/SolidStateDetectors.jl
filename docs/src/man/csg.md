# Constructive Solid Geometry (CSG)

All objects are defined through Constructive Solid Geometry (CSG),
where complex geometries can be constructed by combining simple volume primitives (e.g. `Tube`) through [Boolean operators](@ref) and transformed using [Transformations](@ref). 

The primitives which can be used are shown under [Volume Primitives](@ref) together with how they
can be specified in the configuration files.


## Volume Primitives

````@example primitives
using SolidStateDetectors
import SolidStateDetectors.ConstructiveSolidGeometry as CSG
using Plots
T = Float64;
nothing #hide
````

### List of YAML example configuration files for Primitives
Under `SolidStateDetectors.jl/examples/example_primitive_files` there
are some examples how to define the different primitives
via the YAML format:

````@example primitives
path_to_example_primitives_config_files = joinpath(dirname(dirname(pathof(SolidStateDetectors))), "examples", "example_primitive_files")
example_primitives_config_filenames = readdir(path_to_example_primitives_config_files)
for fn in example_primitives_config_filenames
    println(fn)
end
````

### Box

````@example primitives
cfn = joinpath(path_to_example_primitives_config_files, "Box.yaml")
print(open(f -> read(f, String), cfn))
````

Load the primitive from the configuration file via `CSG.Geometry`

````@example primitives
box = CSG.Geometry(T, cfn)
plot(box)
````

### Cone
#### Tube

````@example primitives
cfn = joinpath(path_to_example_primitives_config_files, "Cone_tube.yaml")
print(open(f -> read(f, String), cfn))
````

Load the primitive from the configuration file via `CSG.Geometry`

````@example primitives
cone = CSG.Geometry(T, cfn)
plot(cone)
````

#### VaryingTube

````@example primitives
cfn = joinpath(path_to_example_primitives_config_files, "Cone.yaml")
print(open(f -> read(f, String), cfn))
````

Load the primitive from the configuration file via `CSG.Geometry`

````@example primitives
cone = CSG.Geometry(T, cfn)
plot(cone)
````

#### Polycone


````@example primitives
cfn = joinpath(path_to_example_primitives_config_files, "Polycone.yaml")
print(open(f -> read(f, String), cfn))
````

Load the primitive from the configuration file via `CSG.Geometry`

````@example primitives
polycone = CSG.Geometry(T, cfn)
plot(polycone)
````


### Ellipsoid
#### Sphere

````@example primitives
cfn = joinpath(path_to_example_primitives_config_files, "Ellipsoid_full_sphere.yaml")
print(open(f -> read(f, String), cfn))
````

Load the primitive from the configuration file via `CSG.Geometry`

````@example primitives
ellipsoid = CSG.Geometry(T, cfn)
plot(ellipsoid)
````

### Torus

````@example primitives
cfn = joinpath(path_to_example_primitives_config_files, "Torus.yaml")
print(open(f -> read(f, String), cfn))
````

Load the primitive from the configuration file via `CSG.Geometry`

````@example primitives
torus = CSG.Geometry(T, cfn)
plot(torus, zlims = [-6,6], camera = (40, 55))
````

### Prism
#### Hexagonal Prism

````@example primitives
cfn = joinpath(path_to_example_primitives_config_files, "RegularPrism_hexagon.yaml")
print(open(f -> read(f, String), cfn))
````

Load the primitive from the configuration file via `CSG.Geometry`

````@example primitives
prism = CSG.Geometry(T, cfn)
plot(prism)
````


## Boolean operators

The Boolean operators are `union`, `difference` and `intersection`:

### Union

A `union` of two objects `A` and `B` is defined as the set of points that are in at least one of either `A` or `B`.
In the configuration files, it is defined using the `union` field, followed by an array of entries to construct the union, e.g.
```yaml
union: # A || B
  - tube: # A
      r: 2
      h: 1
  - tube: # B
      r: 1
      h: 1.5
      origin: 
        z: 0.5
```
![CSGUnion](../assets/CSGUnion.png)

If more than two geometries are passed, the `union` is constructed from all of them.


### Difference


A `difference` of two objects `A` and `B` is defined as the set of points that are in `A` but not in `B`. Note that `B` is treated as open primitive. This means that points which are in `A` and on the surface of `B` will still be in the `difference` of `A` and `B`.
In the configuration files, it is defined using the `difference` field, followed by an array of entries. The first entry of the array is the main geometry, from which all following geometry entries are subtracted, e.g.
```yaml
difference: # A && !B
  - tube: # A
      r: 2
      h: 1
  - tube: # B
      r: 1
      h: 1.1
```
![CSGDifference](../assets/CSGDifference.png)

Keep in mind that to discard the part of the surface of `A` which is on the surface of `B`, `B` should be chosen slightly bigger than `A`.

If more than two geometries are passed, all entries starting from the second will be subtracted from the first.


### Intersection

An `intersection` of two objects `A` and `B` is defined as the set of points that are both in `A` and in `B`.
In the configuration files, it is defined using the `intersection` field, followed by an array of entries to construct the intersection, e.g.
```yaml
intersection: # A && B
  - tube: # A
      r: 2
      h: 1
  - tube: # B
      r: 1
      h: 1.5
      origin: 
        z: 0.5
```
![CSGIntersection](../assets/CSGIntersection.png)

If more than two geometries are passed, the `intersection` is constructed from all of them.


## Transformations

All [Volume Primitives](@ref) are defined such that they are centered around the origin of the coordinate system. They can be rotated in their local coordinate system and translated to their final position in the global coordinate system.

There are two possibilities two rotate and translate volume primitives. One is to define the rotation/translation inside the primitive definition. The other one is to define it for complete sets like detectors, unions, etc.

### Rotations

Rotations are defined in the configuration files by either a 3$\times$3 rotation matrix or a set of angles with respective rotation axes are required.

If a `tube` is to be rotated 45° around the `x` axis, the rotation is parsed to the primitive as additional `rotation` field in the primitive definition.
```yaml
tube:
  r: 1
  h: 1
  rotation:
    X: 45°
```
In this case, the rotation around the `x` axis by 45° is parsed in the format seen above: the name of the axis to be rotated around is given as (upper-case) letter, followed by the rotation angle. If no units are given to the rotation angle, it will be parsed with `units.angle`.

If the rotation is to be described as multiple subsequent rotations, it can be passed using an array of angles with respective axis description as field name, e.g.
```yaml
tube:
  r: 1
  h: 1
  rotation:
    XZ: [45°, 30°]
```
will first rotate the tube 45° around the `x` axis, followed by a 30° rotation around the `z` axis.

Alternatively, a full 3$\times$3 matrix can be passed using the `M` field, e.g.
```yaml
tube:
  r: 1
  h: 1
  rotation:
    M: [1, 0, 0, 0, 0, -1, 0, 1, 0]
```
will transform the primitive using the rotation matrix
```math
\left[ \begin{array}{ccc}1&0&0\\0&0&-1\\0&1&0\end{array}\right]
```
which corresponds to `X: 90°`.

Alternative naming for the `rotation` field can be `rotate`.

!!! note
    If $\varphi$ is specified for a certain interval, the interval is internally converted into a rotation, $R_{\varphi}$, of the primitive.
    This rotation is applied prior to the rotation specified in the rotation field, $R$. 
    Thus, internally, both rotations are calculated into the final rotation matrix of the primitive, $R_f$, via $R_f = R \cdot R_{\varphi}$. 


### Translations

Translations are defined in the configuration files through a Cartesian vector.

If a `tube` is translated 1cm along the `x` axis, the translation is parsed to the primitive as additional `origin` field in the primitive definition.
```yaml
tube:
  r: 1
  h: 1
  origin:
    x: 1cm
```
The Cartesian vector can also be passed as a vector, i.e.
```yaml
tube:
  r: 1
  h: 1
  origin: [1cm, 0, 0]
```
If no units are given, the translation is parsed in units of `units.length`.

Alternative naming for the `origin` field can be `translate` and `translation`.



### Combination of Transformations

If both a rotation and translation are defined in a primitive definition, it is first rotated and then translated.
```yaml 
tube:
  r: 1
  h: 1
  rotation:
    X: 45°
  origin: 
    z: 1cm
```
would first rotate the `tube` by 45° around the `x` axis before translating it 1cm along the `z` axis.
  

### Transforming Sets

If a `union` of two primitives should be transformed all together, the transformation can also be defined by nesting the `union` into a `translate` with respective information, e.g.
```yaml
translate: 
  z: 1
  union:
    - tube: 
        r: 1
        h: 1
    - box: 
        widths: [1,1,1]
```

Same applies for rotations or other sets, e.g.
```yaml
rotate: 
  X: 45°
  difference:
    - tube: 
        r: 1
        h: 1
    - box: 
        widths: [1,1,1]
```