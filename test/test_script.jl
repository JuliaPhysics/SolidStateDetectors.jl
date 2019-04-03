#=
    This is script to test the hole simulation chain with all standard detector types. 
    The script also produces some output plots but only very basic ones, since this script
    is only to test the core functionality of the package.
=#

outputdir = joinpath(ENV["HOME"], "tmp/solidstatedetectors.jl/")
mkpath(outputdir)
@info "Test output dir: $outputdir"

@info "Loading packages"
using Plots; pyplot()
using SolidStateDetectors

ivc = SolidStateDetector(SSD_examples[:InvertedCoax])
coax = SolidStateDetector(SSD_examples[:Coax])
bege = SolidStateDetector(SSD_examples[:BEGe])
cgd = SolidStateDetector(SSD_examples[:CGD])

plot()
for key in keys(SSD_examples)
    @info "Now test detector type: $key"
    det = SolidStateDetector(SSD_examples[key])

    g = Grid(det)

    Init_setup = SolidStateDetectors.PotentialSimulationSetup(det);
    plot(Init_setup, size = (1200, 1200))
    savefig(joinpath(outputdir, "$(key)_init_setup"))
    
    for nrefs in [0, 1, 2, 3]
        setup = calculate_electric_potential(det, max_refinements = nrefs)
        plot(setup, size = (1200, 1200))
        savefig(joinpath(outputdir, "$(key)_setup_$(nrefs)_refinements"))
    end    
end

@info "Finished testing."
@info "Test output saved in: $outputdir"
