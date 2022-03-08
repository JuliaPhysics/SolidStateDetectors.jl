# Use
#
# DOCUMENTER_DEBUG=true GKSwstype=100 julia --color=yes make.jl local [linkcheck] [fixdoctests]
#
# for local builds.

using Pkg
using Documenter
using Literate
using Plots
using SolidStateDetectors
using SolidStateDetectors.ConstructiveSolidGeometry

function fix_literate_output(content)
    content = replace(content, "EditURL = \"@__REPO_ROOT_URL__/\"" => "")
    return content
end

gen_content_dir = joinpath(@__DIR__, "src")
tutorial_src = joinpath(@__DIR__, "src", "tutorial_lit.jl")
Literate.markdown(tutorial_src, gen_content_dir, name = "tutorial", documenter = true, credit = true, postprocess = fix_literate_output)
#Literate.markdown(tutorial_src, gen_content_dir, name = "tutorial", codefence = "```@repl tutorial" => "```", documenter = true, credit = true)
Literate.notebook(tutorial_src, gen_content_dir, execute = false, name = "ssd_tutorial", documenter = true, credit = true)
Literate.script(tutorial_src, gen_content_dir, keep_comments = false, name = "ssd_tutorial", documenter = true, credit = false)

makedocs(
    sitename = "SolidStateDetectors.jl",
    modules = [SolidStateDetectors, SolidStateDetectors.ConstructiveSolidGeometry],
    pages = [
        "Home" => "index.md",
        "Manual" => Any[
            "Installation" => "man/installation.md",
            "Configuration Files" => "man/config_files.md",
            "Constructive Solid Geometry (CSG)" => "man/csg.md",
            "Grids" => "man/Grids.md",
            "Electric Potential" => "man/electric_potential.md",
            "Electric Field" => "man/electric_field.md",
            "Charge Drift" => "man/charge_drift.md",
            "Weighting Potentials" => "man/weighting_potentials.md",
            "Capacitances" => "man/capacitances.md",
            "IO" => "man/IO.md",
            "Plotting" => "man/plotting.md",
        ],
        "Tutorial" => "tutorial.md",
        "API" => "api.md",
        "LICENSE" => "LICENSE.md",
    ],
    doctest = ("fixdoctests" in ARGS) ? :fix : true,
    format = Documenter.HTML(canonical = "https://JuliaPhysics.github.io/SolidStateDetectors.jl/stable/", prettyurls = !("local" in ARGS)),
    linkcheck = ("linkcheck" in ARGS),
    strict = !("local" in ARGS),
)

deploydocs(
    repo = "github.com/JuliaPhysics/SolidStateDetectors.jl.git",
    devbranch = "main",
    devurl = "main",
    forcepush = true,
    push_preview = true,
)
