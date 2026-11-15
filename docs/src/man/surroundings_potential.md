# Potential of the Surroundings

A simulation only covers a finite volume, the world (see [Grid](@ref)). Everything outside it, e.g. the cryostat around the detector, is represented by a single voltage: the **potential of the surroundings**, `0 V` (grounded) by default.

Whether it matters depends on the boundary conditions of the world (see [Grid Boundary Conditions](@ref)):
- `fixed`: the edge of the world behaves like a metal wall held at this potential.
- `infinite`: the potential gradually approaches this potential beyond the edge of the world.
- `reflecting`, `periodic`: the surroundings have no effect.

## When to set it

**Almost never.** By default, the surroundings are at `0 V` (grounded), like in most experiments.

| The environment (e.g. cryostat) is... | surrounding potential |
|:---|:---|
| grounded (`0 V`) | `0 V`, the default: nothing to set |
| held at another voltage, e.g. `+500 V` | `calculate_electric_potential!(sim, surroundings_potential = 500u"V")` |

This does **not** depend on the voltages of the contacts. For example, with the contacts at `-1500 V` and `+1500 V` and a grounded cryostat:
```julia
sim = Simulation(SSD_examples[:InvertedCoax])
sim.detector = SolidStateDetector(sim.detector, contact_id = 1, contact_potential = -1500)
sim.detector = SolidStateDetector(sim.detector, contact_id = 2, contact_potential = +1500)
calculate_electric_potential!(sim)   # surroundings at 0 V, the default
```

Once set, the value is remembered by `sim` and used by all later calculations of its electric potential.
It also works with [`simulate!`](@ref).


## Only voltage differences matter

A voltage is always measured with respect to something. If you describe the same setup with a different "zero",
you have to change **every** voltage by the same amount, **including the surroundings**. Otherwise, a different setup is described.

Example: the Inverted Coax detector with its point contact grounded and its mantle at `2016 V`:

| | Point contact | Mantle | Surroundings | Same setup as A? |
|:---|:---:|:---:|:---:|:---|
| **A** (original) | `0 V` | `2016 V` | `0 V` | – |
| **B** (all voltages `-2016 V`) | `-2016 V` | `0 V` | `-2016 V` | **Yes**: same results, the potential is just `2016 V` lower everywhere |
| **C** (only the contacts `-2016 V`) | `-2016 V` | `0 V` | `0 V` | **No**: now the mantle, not the point contact, has the same voltage as the surroundings |

This can be checked directly:
```julia
using SolidStateDetectors

function is_depleted_ivc(boundary, V_point, V_mantle, V_surroundings)
    cfg = SolidStateDetectors.parse_config_file(SSD_examples[:InvertedCoax])
    cfg["grid"]["axes"]["r"]["boundaries"] = boundary
    cfg["grid"]["axes"]["z"]["boundaries"] = boundary
    sim = Simulation{Float64}(cfg)
    sim.detector = SolidStateDetector(sim.detector, contact_id = 1, contact_potential = V_point)
    sim.detector = SolidStateDetector(sim.detector, contact_id = 2, contact_potential = V_mantle)
    calculate_electric_potential!(sim, depletion_handling = true, surroundings_potential = V_surroundings, verbose = false)
    is_depleted(sim.point_types)
end

for boundary in ("fixed", "reflecting", "inf")
    A = is_depleted_ivc(boundary,     0, 2016,     0)
    B = is_depleted_ivc(boundary, -2016,    0, -2016)
    C = is_depleted_ivc(boundary, -2016,    0,     0)
    println(boundary, ":  A = $A,  B = $B,  C = $C")
end
```

| Boundary | A (original) | B (all voltages shifted) | C (only contacts shifted) |
|:---|:---:|:---:|:---:|
| `fixed` | depleted | depleted | **not depleted** |
| `reflecting` | depleted | depleted | depleted |
| `inf` | depleted | depleted | **not depleted** |

- **B always gives the same result as A**: it is the same setup, written with a different zero.
- **C can differ** for `fixed` and `inf`, because it is a different setup. At `2016 V`, this detector is just barely depleted, so this small change matters.
  For `reflecting`, the surroundings play no role, so C gives the same result as A.


## How much the surroundings matter

The closer the edge of the world is to the detector, the more the surroundings affect the electric potential inside it.
This is strongest for `infinite` boundaries, which only approximate open space at a finite distance.

A wrong value can therefore change the results, e.g. the depletion voltage.
For example, a detector with a grounded mantle and its point contact at `-3000 V` is described correctly with the default (`0 V`).
Setting `surroundings_potential = -3000` instead would put the cryostat at the voltage of the point contact, which is a different setup.

A larger world reduces this sensitivity, but does not replace the correct value: use the potential your environment really has.


## Notes

- **Weighting potentials** always use grounded surroundings. `surroundings_potential` has no effect on them.
- **Saved simulations** do not store `surroundings_potential`. After loading one, pass it again before recalculating the electric potential.
- **Step by step**: if you use [`apply_initial_state!`](@ref) and `SolidStateDetectors.update_till_convergence!` directly,
  set the value first with `sim.world = SolidStateDetectors.World(sim.world, surroundings_potential = 500)`.
- **Depletion voltage**: if no contact is at `0 V`, pass `contact_id` and the voltage range to [`estimate_depletion_voltage`](@ref).
