using Plots; pyplot()
using SolidStateDetectors

ivc = SolidStateDetector(SSD_examples[:InvertedCoax])
coax = SolidStateDetector(SSD_examples[:Coax])
bege = SolidStateDetector(SSD_examples[:BEGe])
cgd = SolidStateDetector(SSD_examples[:CGD])

plot()
for key in keys(SSD_examples)
    @info "Now: $key"
    det = SolidStateDetector(SSD_examples[key])

    g = Grid(det)

    Init_setup = SolidStateDetectors.PotentialSimulationSetup(det);
    plot(Init_setup, size = (1200, 1200))
    savefig("/home/iwsatlas1/lhauert/tmp/$(key)_init_setup")
    
    for nrefs in [0, 1, 2, 3]
        setup = calculate_electric_potential(det, max_refinements = nrefs)
        plot(setup, size = (1200, 1200))
        savefig("/home/iwsatlas1/lhauert/tmp/$(key)_setup_$(nrefs)_refinements")
    end    
end