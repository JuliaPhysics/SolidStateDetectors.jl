# This file is a part of SolidStateDetectors.jl, licensed under the MIT License (MIT).


const elementary_charge = Float64(1.602176487e-19)
const ϵ0  = Float64(8.8541878176e-12)

const material_properties = Dict{Symbol, NamedTuple}()

abstract type AbstractDriftMaterial end

material_properties[:Vacuum] = (
    E_ionisation = 0.0u"eV",
    f_fano = 0.0,
    ϵ_r = 1.0,
    ρ = 0u"g*cm^-3",
    name = "Vacuum",
    ml = 1.0,
    mt = 1.0
)

abstract type HPGe <: AbstractDriftMaterial end 
Symbol(::Type{HPGe}) = :HPGe
material_properties[:HPGe] = (
    E_ionisation = 2.95u"eV",
    f_fano = 0.129, # https://doi.org/10.1103/PhysRev.163.238
    ϵ_r = 16.0,
    ρ = 5.323u"g*cm^-3",
    name = "High Purity Germanium",
    ml = 1.64,
    mt = 0.0819,
    diffusion_fieldvector_electrons = 700,
    diffusion_fieldvector_holes = 250
)


# These values might just be approximations
abstract type Si <: AbstractDriftMaterial end
Symbol(::Type{Si}) = :Si
material_properties[:Si] = (
    E_ionisation = 3.62u"eV",
    f_fano = 0.11,
    ϵ_r = 11.7,
    ρ = 2.3290u"g*cm^-3",
    name = "Silicon",
    mt = 0.98,
    ml = 0.19
)

material_properties[:Al] = (
    name = "Aluminium",
    ϵ_r = 10.8, # Aluminium Foil
    ρ = 2.6989u"g*cm^-3"
)

material_properties[:LAr] = (
    name = "liquid Argon",
    E_ionisation = 0u"eV",
    f_fano = 0.107,
    ϵ_r = 1.505,
    ρ = 1.396u"g*cm^-3"
)

material_properties[:Cu] = (
    name = "Copper",
    ϵ_r = 20,
    ρ = 8.96u"g*cm^-3"
)

material_properties[:CdZnTe] = (
    name = "Cadmium zinc telluride",
    E_ionisation = 4.64u"eV",
    f_fano = 0.089, # https://doi.org/10.1557/PROC-487-101
    ϵ_r = 10.9,
    ρ = 5.78u"g*cm^-3"
)
# Add new materials above this line
# and just put different spellings into the dict `materials` below

materials = Dict{String, Symbol}(
    "HPGe" => :HPGe,
    "vacuum" => :Vacuum,
    "Vacuum" => :Vacuum,
    "Copper" => :Cu,
    "copper" => :Cu,
    "Al"  => :Al,
    "LAr" => :LAr,
    "CZT" => :CdZnTe,
    "Si" => :Si,
    "Co" => begin
        @warn "In v0.9.0, the material 'Co' does not denote copper anymore. Please use 'Cu' for copper."
        :Co
    end
)


for key in keys(material_properties)
    if !haskey(materials, string(key))
        push!(materials, string(key) => key)
    end
end
