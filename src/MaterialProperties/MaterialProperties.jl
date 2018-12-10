# This file is a part of SolidStateDetectors.jl, licensed under the MIT License (MIT).

const material_properties = Dict{Symbol,NamedTuple}()

material_properties[:Vacuum] = (
    E_ionisation = 0.0u"eV",
    f_fano = 0.0,
    ϵ_r = 1.0,
    ρ = 0Unitful.g / (Unitful.cm^3),
    name = "Vacuum"
)

material_properties[:HPGe] = (
    E_ionisation = 2.95u"eV",
    f_fano = 0.13,
    ϵ_r = 16.0,
    ρ = 5.323Unitful.g / (Unitful.cm^3), # u"g/cm^3" throws warnings in precompilation
    name = "High Purity Germanium"
)
