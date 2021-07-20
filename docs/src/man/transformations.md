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
which corresponds to `X: 45°`.

Alternative naming for the `rotation` field can be `rotate`.


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