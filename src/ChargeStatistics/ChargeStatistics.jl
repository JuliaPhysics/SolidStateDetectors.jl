# This file is a part of SolidStateDetectors.jl, licensed under the MIT License (MIT).


"""
    add_fano_noise(E_dep::RealQuantity, E_ionisation::RealQuantity, f_fano::Real)::RealQuantity

Adds Fano noise to an energy deposition `E_dep`, assuming a detector material
ionisation energy `E_ionisation` and a Fano factor `f_fano`.

## Arguments
* `E_dep::RealQuantity`: Energy deposited in a semiconductor material.
* `E_ionisation`: Energy needed to create one electron-hole-pair in the semiconductor material.
* `f_fano`: Fano factor of the material.

## Example 
```julia
add_fano_noise(100u"keV", 2.95u"eV", 0.129)
```

Some material properties are stored in `SolidStateDetectors.material_properties` and can be used here:
```julia 
material = SolidStateDetectors.material_properties[:HPGe]
add_fano_noise(100u"keV", material.E_ionisation, material.f_fano)
```

!!! note
    Using values with units for `E_dep` or `E_ionisation` requires the package [Unitful.jl](https://github.com/PainterQubits/Unitful.jl).
"""
function add_fano_noise end
export add_fano_noise

function add_fano_noise(rng::AbstractRNG, E_dep::RealQuantity, E_ionisation::RealQuantity, f_fano::Real)
    target_unit = unit(E_dep)
    n_expected = uconvert(Unitful.NoUnits, E_dep/E_ionisation)

    # dist = Distributions.Poisson(n_expected)
    # n_observed = n_expected + f_fano * (rand(dist) - n_expected)

    # Use normal distribution instead of Poisson, makes no difference for
    # typical numbers of charges:
    σ = sqrt(f_fano * n_expected)
    n_observed = n_expected + σ * randn(rng)
    n_observed_pos = max(0, n_observed)

    uconvert(target_unit, E_ionisation * n_observed_pos)
end


function add_fano_noise(rng::AbstractRNG, evts::DetectorHitEvents, E_ionisation::RealQuantity, f_fano::Real)
    noisy_edep = deepmap(x -> add_fano_noise(x, E_ionisation, f_fano), evts.edep)

    TypedTables.Table(merge(
        Tables.columns(evts),
        (edep = noisy_edep,)
    ))
end


add_fano_noise(x, E_ionisation::RealQuantity, f_fano::Real) =
    add_fano_noise(Random.GLOBAL_RNG, x, E_ionisation, f_fano)
