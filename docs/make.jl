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

gen_content_dir = joinpath(@__DIR__, "src", "tutorials")
for tut_lit_fn in filter(fn -> endswith(fn, "_lit.jl"), readdir(gen_content_dir))
    lit_src_fn = joinpath(gen_content_dir, tut_lit_fn)
    tut_basename = tut_lit_fn[1:end-7] # remove "_lit.jl"
    Literate.markdown(lit_src_fn, gen_content_dir, name = tut_basename, documenter = true, credit = true, postprocess = fix_literate_output)
    Literate.notebook(lit_src_fn, gen_content_dir, execute = false, name = tut_basename, documenter = true, credit = true)
    Literate.script(lit_src_fn, gen_content_dir, keep_comments = false, name = tut_basename, documenter = true, credit = false)
end


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
        "Tutorials" => [
            "tutorials/complete_simulation_chain_IVC.md",
            "tutorials/custom_impurity_density_pn_junction.md",
        ],
        "API" => "api.md",
        "LICENSE" => "LICENSE.md",
    ],
    doctest = ("fixdoctests" in ARGS) ? :fix : true,
    format = Documenter.HTML(canonical = "https://JuliaPhysics.github.io/SolidStateDetectors.jl/stable/", size_threshold = 1048576, prettyurls = !("local" in ARGS)),
    linkcheck = ("linkcheck" in ARGS),
    warnonly = ("nonstrict" in ARGS),
)

deploydocs(
    repo = "github.com/JuliaPhysics/SolidStateDetectors.jl.git",
    devbranch = "main",
    devurl = "main",
    forcepush = true,
    push_preview = true,
)
