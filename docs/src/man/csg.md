# Constructive Solid Geometry (CSG)

All objects are defined through Constructive Solid Geometry (CSG),
where complex geometries can be constructed by combining simple volume primitives (e.g. `Tube`)
through Boolean operators. 

The primitives which can be used are shown under [Volume Primitives](@ref) together with how they
can be specified in the configuration files.


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


A `difference` of two objects `A` and `B` is defined as the set of points that are in `A` but not in `B`.
In the configuration files, it is defined using the `difference` field, followed by an array of entries. The first entry of the array is the main geometry, from which all following geometry entries are subtracted, e.g.
```yaml
difference: # A && !B
  - tube: # A
      r: 2
      h: 1
  - tube: # B
      r: 1
      h: 1
```
![CSGDifference](../assets/CSGDifference.png)

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