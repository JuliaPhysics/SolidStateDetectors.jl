# Use
#
# DOCUMENTER_DEBUG=true julia --color=yes make.jl local [linkcheck] [fixdoctests]
#
# for local builds.

using Documenter
using Literate
using Plots
pyplot(fmt=:png) 
using SolidStateDetectors

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
    modules = [SolidStateDetectors],
    pages=[
        "Home" => "index.md",
        "Manual" => Any[
            "Installation" => "man/installation.md",
            "Detectors" => Any[
                "Config Files" => "man/config_files.md",
            ],
            "Geometries (CSG)" => Any[
                "CSG" => "man/csg.md",
                "Primitives" => "man/primitives.md",
            ],
            "Electric Potentials" => "man/electric_potentials.md",
            "Weighting Potentials" => "man/weighting_potentials.md",
            "Electric Fields" => "man/electric_fields.md",
            "Drift Fields" => "man/drift_fields.md",
            "IO" => "man/IO.md",
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
    devbranch = "master",
    devurl = "master",
    forcepush = true,
    push_preview = true,
)
