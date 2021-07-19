# This file is a part of SolidStateDetectors.jl, licensed under the MIT License (MIT).


const elementary_charge = Float64(1.602176487e-19)
const ϵ0  = Float64(8.8541878176e-12)

const material_properties = Dict{Symbol, NamedTuple}()

abstract type AbstractDriftMaterial end

material_properties[:Vacuum] = (
    E_ionisation = 0.0u"eV",
    f_fano = 0.0,
    ϵ_r = 1.0,
    ρ = 0Unitful.g / (Unitful.cm^3),
    name = "Vacuum",
    ml = 1.0,
    mt = 1.0
)

abstract type HPGe <: AbstractDriftMaterial end 
Symbol(::Type{HPGe}) = :HPGe
material_properties[:HPGe] = (
    E_ionisation = 2.95u"eV",
    f_fano = 0.13,
    ϵ_r = 16.0,
    ρ = 5.323Unitful.g / (Unitful.cm^3), # u"g/cm^3" throws warnings in precompilation
    name = "High Purity Germanium",
    ml = 1.64,
    mt = 0.0819
)


# These values might just be approximations
abstract type Si <: AbstractDriftMaterial end
Symbol(::Type{Si}) = :Si
material_properties[:Si] = (
    E_ionisation = 3.62u"eV",
    f_fano = 0.11,
    ϵ_r = 11.7,
    ρ = 2.3290Unitful.g / (Unitful.cm^3), # u"g/cm^3" throws warnings in precompilation
    name = "Silicon",
    mt = 0.98,
    ml = 0.19
)

material_properties[:Al] = (
    name = "Aluminium",
    ϵ_r = 10.8, # Aluminium Foil
    ρ = 2.6989Unitful.g / (Unitful.cm^3) # u"g/cm^3" throws warnings in precompilation
)

material_properties[:LAr] = (
    name = "liquid Argon",
    E_ionisation = 0u"eV",
    f_fano = 0.107,
    ϵ_r = 1.505,
    ρ = 1.396Unitful.g / (Unitful.ml) # u"g/cm^3" throws warnings in precompilation
)

material_properties[:Co] = (
    name = "Copper",
    ϵ_r = 20,
    ρ = 8.96Unitful.g / (Unitful.ml) # u"g/cm^3" throws warnings in precompilation
)

material_properties[:CdZnTe] = (
    name = "Cadmium zinc telluride",
    E_ionisation = 4.64u"eV",
    ϵ_r = 10.9,
    ρ = 5.78Unitful.g / (Unitful.ml) # u"g/cm^3" throws warnings in precompilation
)
# Add new materials above this line
# and just put different spellings into the dict `materials` below

materials = Dict{String, Symbol}(
    "HPGe" => :HPGe,
    "vacuum" => :Vacuum,
    "Vacuum" => :Vacuum,
    "Copper" => :Co,
    "copper" => :Co,
    "Al"  => :Al,
    "LAr" => :LAr,
    "CZT" => :CdZnTe,
    "Si" => :Si
)


for key in keys(material_properties)
    if !haskey(materials, string(key))
        push!(materials, string(key) => key)
    end
end
