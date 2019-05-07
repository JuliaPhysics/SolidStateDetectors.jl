# This file is a part of SolidStateDetectors.jl, licensed under the MIT License (MIT).


"""
    add_fano_noise(E_dep::RealQuantity, E_ionisation::RealQuantity, f_fano::Real)::RealQuantity

Add Fano noise to an energy deposition `E_dep`, assuming a detector material
ionisation energy `E_ionisation` and a Fano factor `f_fano`.
"""
function add_fano_noise end
export add_fano_noise

function add_fano_noise(rng::AbstractRNG, E_dep::RealQuantity, E_ionisation::RealQuantity, f_fano::Real)
    target_unit = unit(E_dep)
    n_expected = uconvert(Unitful.NoUnits, E_dep/E_ionisation)

    # dist = Distributions.Poisson(n_expected)
    # n_observed = n_expected + f_fano * (rand(dist) - n_expected)

    # Use normal distribution instead of Poisson, makes no differenc for
    # typical numbers of charges:
    σ = sqrt(n_expected)
    n_observed = n_expected + f_fano * σ * randn(rng)
    n_observed_pos = max(0, n_observed)

    uconvert(target_unit, E_ionisation * n_observed_pos)
end


function add_fano_noise(rng::AbstractRNG, events::DetectorHitEvents, E_ionisation::RealQuantity, f_fano::Real)
    noisy_edep = deepmap(x -> add_fano_noise(x, E_ionisation, f_fano), events.edep)

    TypedTables.Table(merge(
        Tables.columns(events),
        (edep = noisy_edep,)
    ))
end


add_fano_noise(x, E_ionisation::RealQuantity, f_fano::Real) =
    add_fano_noise(Random.GLOBAL_RNG, x, E_ionisation, f_fano)
