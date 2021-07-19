# Constructive Solid Geometry (CSG)

All objects are defined through Constructive Solid Geometry (CSG),
where complex geometries can be build by combining simple volume primitives (e.g. `Box`)
through boolean operators. 

The primitives which can be used are shown under [Volume Primitives](@ref) together with how they
can be specified in the configuration files.
## Boolean operators

The boolean operators are `union`, `difference` and `intersection`:
### Union

![CSGUnion](../assets/CSGUnion.png?raw=true)

```yaml
geometry:
  type: union # A || B
  parts:
    - tube: # A
        r: 0.5
        h: 1
        origin: 
          z: 1
    - tube: # B
        r: 1
        h: 1
```

### Difference

![CSGDifference](../assets/CSGDifference.png?raw=true)

```yaml
geometry: # A && !B
  type: difference
  parts:
    - tube: # A
        r: 0.5
        h: 1.2
    - tube: # B
        r: 1
        h: 1
```

### Intersection

![CSGIntersection](../assets/CSGIntersection.png?raw=true)

```yaml
geometry: # A && B
  type: intersection
  parts:
    - tube: # A
        r: 
          from: 0.5
          to: 1
        h: 1
    - tube: # B
        r: 1
        h: 1
```