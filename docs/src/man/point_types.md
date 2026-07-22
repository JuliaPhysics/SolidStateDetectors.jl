# Point Types

After running `calculate_electric_potential!`, every grid point in `sim.point_types` is assigned a `PointType` value — a `UInt8` bitmask that encodes the physical role of that point in the simulation.
Point types are used throughout the code to decide which points are updated during the successive over-relaxation (SOR, see [Simulation Algorithm](electric_potential.md#Simulation-Algorithm)), whether the detector is depleted, which points belong to the inactive layer, and what volume counts as active.

## Bit Flags

Six single-bit flags are defined, each encoding one physical property:

| Bit | Constant | Hex | Meaning |
|-----|----------|-----|---------|
| 1 | `update_bit` | `0x01` | `1` → point is updated in SOR; `0` → fixed point (contacts, exterior boundaries) |
| 2 | `undepleted_bit` | `0x02` | `1` → point is undepleted; `0` → fully depleted |
| 3 | `pn_junction_bit` | `0x04` | `1` → point lies inside the semiconductor volume |
| 4 | `bulk_bit` | `0x08` | `1` → all points in the 3×3×3 neighbourhood (face-, edge-, and corner-adjacent) also have `pn_junction_bit` and `update_bit` set |
| 5 | `inactive_layer_bit` | `0x10` | `1` → point is part of the inactive (dead) layer |
| 6 | `inactive_contact_bit` | `0x20` | `1` → point is part of the contact adjacent to the inactive layer |

They can be queried using the provided convenience functions or directly with bitwise AND:

```julia
import SolidStateDetectors: update_bit, undepleted_bit, pn_junction_bit,
                             bulk_bit, inactive_layer_bit, inactive_contact_bit,
                             is_fixed_point_type, is_undepleted_point_type,
                             is_pn_junction_point_type, is_in_inactive_layer

pt = sim.point_types.data[i, j, k]   # PointType at one grid point

is_fixed_point_type(pt)               # pt & update_bit         == 0
is_undepleted_point_type(pt)          # pt & undepleted_bit     > 0
is_pn_junction_point_type(pt)         # pt & pn_junction_bit    > 0
is_in_inactive_layer(pt)              # pt & inactive_layer_bit > 0
pt & bulk_bit             > 0         # true if the full 3×3×3 neighbourhood is inside the semiconductor
pt & inactive_contact_bit > 0         # true if the point is part of the doped contact adjacent to the inactive layer
```

## Common Point Type Values

The six flags combine to describe the physical state of each grid point.
The most frequently encountered values are:

| Physical description | Dec | Hex | Active bits |
|----------------------|-----|-----|-------------|
| Contact (fixed) | 0 | `0x00` | None |
| Exterior / vacuum (updateable) | 1 | `0x01` | `update` |
| Semiconductor boundary, depleted | 5 | `0x05` | `pn_junction` + `update` |
| Semiconductor boundary, undepleted | 7 | `0x07` | `pn_junction` + `undepleted` + `update` |
| Active bulk, depleted | 13 | `0x0d` | `bulk` + `pn_junction` + `update` |
| Active bulk, undepleted | 15 | `0x0f` | `bulk` + `pn_junction` + `undepleted` + `update` |
| Inactive layer (boundary) | 23 | `0x17` | `inactive_layer` + `pn_junction` + `undepleted` + `update` |
| Inactive layer | 31 | `0x1f` | `inactive_layer` + `bulk` + `pn_junction` + `undepleted` + `update` |
| Contact adjacent to inactive layer | 32 | `0x20` | `inactive_contact` |


## Information from Point Types

### Depletion and Active Volume

`is_depleted(sim.point_types)` returns `true` when no **bulk** point remains undepleted outside the inactive layer.
Concretely, it checks whether any point simultaneously satisfies all three conditions:

```
bulk_bit = 1    AND    undepleted_bit = 1    AND    inactive_layer_bit = 0
```

If no such point exists, the detector is considered fully depleted.
Points inside the inactive layer (`inactive_layer_bit = 1`) are **excluded** from this check.

The active volume (`get_active_volume(sim.point_types)`) sums the cell volumes of all grid cells satisfying:

```
pn_junction_bit = 1    AND    undepleted_bit = 0    AND    update_bit = 1
```

Because inactive layer cells always carry `undepleted_bit = 1`, they are automatically excluded from the active volume.

### Depletion Voltage Estimation

`estimate_depletion_voltage` applies the same inactive layer exclusion: only points with `pn_junction_bit = 1 AND inactive_layer_bit = 0` (and `bulk_bit = 1` in the final refinement step) are considered.
The returned voltage is therefore the bias at which the **active bulk** outside the inactive layer first becomes fully depleted.

### Inactive Layer Automatic Detection

The inactive layer is the permanently undepleted region that forms near highly-doped contact of a semiconductor detector.

When `calculate_electric_potential!` is run with `depletion_handling = true`, `mark_inactivelayer_bits!` automatically identifies the inactive layer.
Starting from every contact point that carries `inactive_contact_bit`, it scans in both directions along each of the three grid axes.
All consecutive points where `pn_junction_bit`, `undepleted_bit`, and `update_bit` are all set receive the `inactive_layer_bit` flag; the scan along each direction stops as soon as a point fails this condition.

The contact carrying `inactive_contact_bit` is identified through the semiconductor impurity density model (`PtypePNJunctionImpurityDensity`).

### Full Depletion Depth via `get_full_depletion_depth`

`get_full_depletion_depth(sim)` computes the Full Depletion Depth (FDD) — the thickness of the inactive layer — at every point on the inner boundary of the inactive layer.
It can also be called directly on a `PointTypes` object as `get_full_depletion_depth(pt)`.

Two contours are extracted from the point types grid:

- **Inner contour** (`r_inner`, `z_inner`): the inactive-layer / active-bulk transition.
  A grid point belongs to the inner contour when it carries `inactive_layer_bit = 1` **and** has at least one face-adjacent neighbour with `inactive_layer_bit = 0` and `pn_junction_bit = 1` (i.e. a neighbour that is inside the semiconductor but outside the inactive layer).

- **Outer contour** (`r_outer`, `z_outer`): the inner face of the doped contact, where the inactive layer meets the detector boundary.
  Every grid tick inside the doped contact geometry carries `inactive_contact_bit = 1`.
  For each inner contour point the nearest such contact tick (by Euclidean distance in the r-z plane) is selected as `r_outer`/`z_outer`.

In both calculation paths the FDD measures the distance from the inner contour to the inner face of the contact; the contact's own physical thickness is never included:

- **When the model stores a `distance_to_contact` function** (e.g. `ThermalDiffusionLithiumDensity` or `PtypePNJunctionImpurityDensity`): `thickness` is computed by evaluating `distance_to_contact` at the inner contour point. By default this function wraps `ConstructiveSolidGeometry.distance_to_surface` on the contact geometry, giving an exact analytical perpendicular distance to the nearest (inner) surface of the contact regardless of surface orientation and regardless of how thick the contact is.
- **Fallback**: Euclidean distance from the inner contour point to the nearest contact grid tick (r-z plane for cylindrical, 3D for Cartesian). Because the inner contour sits just outside the contact, the nearest tick is on the innermost layer of the contact for standard detector geometries, so the contact thickness is not counted here either.

This perpendicular distance is stored as the `thickness` field of each result entry (see field table below). The Euclidean distance between `(r_inner, z_inner)` and `(r_outer, z_outer)` agrees with `thickness` for axis-aligned surfaces but may differ for diagonal ones:

For **axis-aligned surfaces** (vertical bore wall, horizontal top/bottom face, vertical outer wall) the nearest contact tick lies exactly along the surface normal, so the Euclidean distance from `(r_inner, z_inner)` to `(r_outer, z_outer)` equals `thickness` exactly.
For **non-axis-aligned surfaces** (diagonal, rotated, or curved, e.g. a cone face) the nearest contact tick may not lie along the surface normal, so the Euclidean distance can differ from `thickness` of the FDD (either direction, depending on grid snapping); use the `thickness` field for the correct value.

Each entry in the result is a `NamedTuple` with fields:

| Field | Description |
|-------|-------------|
| `r_inner` / `z_inner` | Coordinates of the inner boundary point (inactive layer / active bulk transition), [m]|
| `r_outer` / `z_outer` | Coordinates of the nearest grid point on the doped contact surface (carries `inactive_contact_bit`), [m]|
| `thickness` | Perpendicular distance from the inner boundary to the contact surface = FDD thickness, [m]|

For Cartesian detectors the fields are `x_inner`, `y_inner`, `z_inner`, `x_outer`, `y_outer`, `z_outer`, `thickness`.

````@setup fdd
using Plots
using SolidStateDetectors
using Unitful
gr(fmt = :png) #hide

T = Float64
sim = Simulation{T}(SSD_examples[:IVCIlayer])
sim.world = typeof(sim.world)(sim.world.intervals, ntuple(_ -> T(9e-5), 3))  # surface refinement spacing
calculate_electric_potential!(sim, depletion_handling = true, use_nthreads = 4)
all_fdd = get_full_depletion_depth(sim)                             # full contour, all coordinates
fdd_at_r = get_full_depletion_depth(sim; r = 5u"mm")                # only inner-contour points nearest r = 5 mm
````

````@example fdd
t_mm = [pt.thickness * 1000 for pt in all_fdd]

p1 = plot(sim.point_types, all_point_types = true,
          title = "FDD inner contour (red)", legend = false, colorbar = true)
scatter!(p1,
    [pt.r_inner for pt in all_fdd], [pt.z_inner for pt in all_fdd],
    markershape = :circle, markersize = 1,
    color = :red, markerstrokewidth = 0)

p2 = scatter(
    [pt.r_inner * 1000 for pt in all_fdd],
    [pt.z_inner * 1000 for pt in all_fdd],
    marker_z = t_mm,
    color = :heat, colorbar = true, colorbar_title = "thickness [mm]",
    clims = extrema(t_mm),
    xlabel = "r [mm]", ylabel = "z [mm]",
    title = "FDD thickness",
    markersize = 3, markerstrokewidth = 0,
    xlims = (0, 40), legend = false, aspect_ratio = :equal
)

plot(p1, p2, layout = (1, 2), size = (900, 450))
````