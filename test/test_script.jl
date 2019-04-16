#=
    This is a script to test the hole simulation chain with all standard detector types. 
    The script also produces some output plots but only very basic ones, since it
    is only to test the core functionality of the package. (So not detector specific plots)
=#

outputdir = joinpath(ENV["HOME"], "tmp/test_solidstatedetectors.jl/")
mkpath(outputdir)
@info "Test output dir: $outputdir"

@info "Loading packages"
using Plots; pyplot()
using SolidStateDetectors

ivc = SolidStateDetector(SSD_examples[:InvertedCoax])
coax = SolidStateDetector(SSD_examples[:Coax])
bege = SolidStateDetector(SSD_examples[:BEGe])
cgd = SolidStateDetector(SSD_examples[:CGD])

T = Float32
@info "Testing now for Float32:"

plot() # creates a plot so that the plots during the following loop pop up.

for key in [:InvertedCoax, :Coax, :BEGe, :CGD]
# for key in keys(SSD_examples)
    @info "Now test detector type: $key"
    
    det = SolidStateDetector{T}(SSD_examples[key])
    S = SSD.get_coordinate_system(det)
    
    setup = SSDSetup(det);
   
    SSD.apply_initial_state!(setup)
    plot(setup.electric_potential)
    savefig(joinpath(outputdir, "$(key)_0_init_setup"))

    for nrefs in [0, 1, 2]
        SSD.calculate_electric_potential!(setup, max_refinements = nrefs)
        plot(setup.electric_potential, size = (1200, 1200))
        savefig(joinpath(outputdir, "$(key)_1_setup_$(nrefs)_refinements"))
    end    
    SSD.calculate_electric_potential!(setup, max_refinements = 3)
    plot(setup.electric_potential)
    savefig(joinpath(outputdir, "$(key)_1_setup_$(3)_refinements"))

    n_contacts = length(setup.detector.contacts)
    for contact in setup.detector.contacts
        SSD.calculate_weighting_potential!(setup, contact.id, max_refinements = 1)
    end

    # plot( # does not work for :Cartesian yet
    #     plot(setup.weighting_potentials[1]),
    #     plot(setup.weighting_potentials[2]),
    #     size = (1000, 1000), layout = (1, 2)
    # )

    SSD.calculate_electric_field!(setup)

    plot( setup.electric_field.grid[1], setup.electric_field.grid[3], SSD.get_electric_field_strength(setup.electric_field)[:, div(length(setup.electric_field.grid[2].ticks), 2), :]', 
          st=:heatmap, clims = (0, 200000), title = "Electric Field Streng [V / m]", xlabel = "x / m", ylabel = "x / m", aspect_ratio = 1, size = (900, 900))
    savefig(joinpath(outputdir, "$(key)_2_Electric_Field_strength"))

    SSD.set_charge_drift_model!(setup, ADLChargeDriftModel())

    SSD.apply_charge_drift_model!(setup)

    pos = if S == :Cartesian
        CartesianPoint{T}[ CartesianPoint{T}( 0.006, 0.00, 0.005  ) ] # this point should be inside all test detectors
    elseif S == :Cylindrical
        CylindricalPoint{T}[ CylindricalPoint{T}( 0.006, 0.0, 0.005  ) ] # this point should be inside all test detectors
    end
    @assert in(pos[1], setup.detector) "Test point $(pos[1]) not inside the detector $(key)."
    begin
        drift_paths = SSD.drift_charges(setup, CartesianPoint.(pos));
        plot(setup.detector)
        plot!(drift_paths)
    end
    savefig(joinpath(outputdir, "$(key)_3_charge_drift"))
end

@info "Finished testing."
@info "Test output saved in: $outputdir"
