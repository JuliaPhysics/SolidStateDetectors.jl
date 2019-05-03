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

T = Float32
@info "Testing now for Float32:"

plot() # creates a plot so that the plots during the following loop pop up.

key = :InvertedCoax

for key in  [:InvertedCoax, :BEGe, :Coax, :CGD]
# for key in keys(SSD_examples)
    @info "Now test detector type: $key"

    det = SolidStateDetector{T}(SSD_examples[key])
    S = SSD.get_coordinate_system(det)

    simulation = Simulation(det);

    SSD.apply_initial_state!(simulation)
    p = if key == :CGD
        plot(simulation.electric_potential, y = 0.002)
    else
        plot(simulation.electric_potential)
    end
    savefig(joinpath(outputdir, "$(key)_0_init_setup"))

    for nrefs in [0, 1, 2]
        SSD.calculate_electric_potential!(simulation, max_refinements = nrefs)
        p = if key == :CGD
            plot(simulation.electric_potential, y = 0.002)
        else
            plot(simulation.electric_potential)
        end
        savefig(joinpath(outputdir, "$(key)_1_Electric_Potential_$(nrefs)_refinements"))
    end
    SSD.calculate_electric_potential!(simulation, max_refinements = 3)
    p = if key == :CGD
        plot(simulation.electric_potential, y = 0.002)
    else
        plot(simulation.electric_potential)
    end
    savefig(joinpath(outputdir, "$(key)_1_Electric_Potential_$(3)_refinements"))

    n_contacts = length(simulation.detector.contacts)
    for contact in simulation.detector.contacts
        SSD.calculate_weighting_potential!(simulation, contact.id, max_refinements = key == :Coax ? 1 : 2)
    end
    wp_plots = if key != :CGD
        [ plot(simulation.weighting_potentials[contact.id]) for contact in simulation.detector.contacts ]
    else
        [ plot(simulation.weighting_potentials[contact.id], y = 0.002) for contact in simulation.detector.contacts ]
    end
    plot( wp_plots..., size = (1000, 1000))
    savefig(joinpath(outputdir, "$(key)_2_Weighting_Potentials"))

    SSD.calculate_electric_field!(simulation)

    plot( simulation.electric_field.grid[1], simulation.electric_field.grid[3], SSD.get_electric_field_strength(simulation.electric_field)[:, div(length(simulation.electric_field.grid[2].ticks), 2), :]',
          st=:heatmap, title = "Electric Field Streng [V / m]", xlabel = "x / m", ylabel = "x / m", aspect_ratio = 1, size = (900, 900))
    savefig(joinpath(outputdir, "$(key)_3_Electric_Field_strength"))

    if S == :cylindrical
        plot_electric_field(simulation, φ=deg2rad(30), spacing = 8)
        savefig(joinpath(outputdir, "$(key)_3_1_Electric_Field_Lines"))
    end

    SSD.set_charge_drift_model!(simulation, ADLChargeDriftModel())

    SSD.apply_charge_drift_model!(simulation)

    pos = if key == :InvertedCoax
        CylindricalPoint{T}[ CylindricalPoint{T}( 0.02, deg2rad(10), 0.025 ) ]
    elseif key == :CGD
        CartesianPoint{T}[ CartesianPoint{T}( 0.006, 0.005, 0.005  ) ] # this point should be inside all test detectors
    elseif key == :BEGe
        CylindricalPoint{T}[ CylindricalPoint{T}( 0.016, deg2rad(10), 0.015  ) ] # this point should be inside all test detectors
    elseif key == :Coax
        CylindricalPoint{T}[ CylindricalPoint{T}( 0.016, deg2rad(10), 0.005  ) ] # this point should be inside all test detectors
    end
    energy_depos = T[1460]
    @assert in(pos[1], simulation.detector) "Test point $(pos[1]) not inside the detector $(key)."

    event = SSD.Event(pos, energy_depos)
    SSD.simulate!(event, simulation)

    plot(simulation.detector)
    plot!(event.drift_paths)
    savefig(joinpath(outputdir, "$(key)_4_charge_drift"))

    # signals[:, 2] *= -1
    plot(event.signals, size = (1200, 600), lw = 1.5)
    savefig(joinpath(outputdir, "$(key)_5_induced_signals"))

end

@info "Finished testing."
@info "Test output saved in: $outputdir"
