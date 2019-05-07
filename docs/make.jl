# Use
#
# DOCUMENTER_DEBUG=true julia --color=yes make.jl local [fixdoctests]
#
# for local builds.

using Documenter
using Pkg
using Plots
pyplot(fmt=:svg)
using SolidStateDetectors

makedocs(
    sitename = "SolidStateDetectors.jl",
    modules = [SolidStateDetectors],
    pages=[
        "Home" => "index.md",
        "Manual" => Any[
            "Installation" => "man/installation.md",
            "Detectors" => "man/detector_geometries.md",
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
    format = Documenter.HTML(canonical = "https://JuliaHEP.github.io/SolidStateDetectors.jl/stable/", prettyurls = !("local" in ARGS))
)

deploydocs(
    repo = "github.com/JuliaHEP/SolidStateDetectors.jl.git",
    devbranch = "master",
    devurl = "master",
)
