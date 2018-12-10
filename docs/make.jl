# Use
#
#     DOCUMENTER_DEBUG=true julia --color=yes make.jl local [fixdoctests]
#
# for local builds.

using Documenter
using SolidStateDetectors

makedocs(
    sitename = "SolidStateDetectors.jl",
    modules = [SolidStateDetectors],
    format = :html,
    pages=[
        "Home" => "index.md",
        "Manual" => Any[
            "Installation" => "man/installation.md",
            "Detectors" => "man/detectors.md",
            "Electric Potentials" => "man/electric_potentials.md",
            "Weighting Potentials" => "man/weighting_potentials.md",
            "Electric Fields" => "man/electric_fields.md",
            "Drift Fields" => "man/drift_fields.md",
            "IO" => "man/IO.md",
        ],
        "API" => "api.md",
        "LICENSE" => "LICENSE.md",
    ],
    doctest = ("fixdoctests" in ARGS) ? :fix : true,
    html_prettyurls = !("local" in ARGS),
    html_canonical = "https://JuliaHEP.github.io/SolidStateDetectors.jl/stable/",
)

deploydocs(
    repo = "github.com/JuliaHEP/SolidStateDetectors.jl.git"
)
