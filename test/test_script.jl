#=
    This is a script to test the whole simulation chain with all standard detector types.
    The script also produces some output plots but only very basic ones, since it
    is only to test the core functionality of the package. (So not detector specific plots)
=#

outputdir = joinpath(ENV["HOME"], "tmp/test_solidstatedetectors.jl/")
mkpath(outputdir)
@info "Test output dir: $outputdir"

@info "Loading packages"
using Plots
using SolidStateDetectors; SSD = SolidStateDetectors
using LinearAlgebra

T = Float32
@info "Testing now for Float32:"

plot() # creates a plot so that the plots during the following loop pop up.

key = :CGD

for key in  [:InvertedCoax, :BEGe, :Coax, :CGD, :Spherical]

    @info "Now test detector type: $key"

    sim = Simulation(SSD_examples[key]);
    S = SSD.get_coordinate_system(sim)

    apply_initial_state!(sim, ElectricPotential)
    p = if S == SSD.Cartesian
        plot(
            plot(sim.electric_potential, y = 0.002),
            plot(sim.q_eff_imp, y = 0.002),
            plot(sim.point_types, y = 0.002),
            size = (1000, 600), layout= (1, 3)
        )
    else
        plot(
            plot(sim.electric_potential),
            plot(sim.q_eff_imp),
            plot(sim.point_types),
            size = (1000, 600), layout= (1, 3)
        )
    end
    savefig(joinpath(outputdir, "$(key)_0_init_setup"))

    nrefs = if key == :InvertedCoax
        0:3
    elseif key == :Spherical
        0:3
    elseif key == :CGD
        0:3
    else
        0:1
    end
    for nref in nrefs
        SSD.update_till_convergence!(sim, ElectricPotential)
        p = if S == SSD.Cartesian
            plot(
                plot(sim.electric_potential, y = 0.002),
                plot(sim.q_eff_imp, y = 0.002),
                plot(sim.point_types, y = 0.002),
                size = (1000, 600), layout= (1, 3)
            )
        else
            plot(
                plot(sim.electric_potential),
                plot(sim.q_eff_imp),
                plot(sim.point_types),
                size = (1000, 600), layout= (1, 3)
            )
        end
        savefig(joinpath(outputdir, "$(key)_1_Electric_Potential_$(nref)_refinements"))
        if nref != nrefs[end]
            SSD.refine!(sim, ElectricPotential, (100, 100, 100), (1e-4, 1e-4, 1e-4), update_other_fields = true)
        end
        @show size(sim.electric_potential.grid)
    end

    for contact in sim.detector.contacts
        calculate_weighting_potential!(sim, contact.id, refinement_limits = key == :Coax ? missing : [0.2], verbose = true)
    end
    wp_plots = if S != SSD.Cartesian
        [ plot(sim.weighting_potentials[contact.id]) for contact in sim.detector.contacts ]
    else
        [ plot(sim.weighting_potentials[contact.id], y = 0.002) for contact in sim.detector.contacts ]
    end
    plot( wp_plots..., size = (1000, 1000))
    savefig(joinpath(outputdir, "$(key)_2_Weighting_Potentials"))

    calculate_electric_field!(sim)

    plot( sim.electric_field.grid[1], sim.electric_field.grid[3], norm.(sim.electric_field)[:, div(length(sim.electric_field.grid[2].ticks), 2), :]',
          st=:heatmap, title = "Electric Field Streng [V / m]", xlabel = "x / m", ylabel = "x / m", aspect_ratio = 1, size = (900, 900))
    savefig(joinpath(outputdir, "$(key)_3_Electric_Field_strength"))

    if S == SSD.Cylindrical
        plot(sim.electric_field, φ = 0, spacing = 3.0)
        savefig(joinpath(outputdir, "$(key)_3_1_Electric_Field_Lines"))
    else
        plot(sim.electric_field, y = 0, spacing = 3.0)
        savefig(joinpath(outputdir, "$(key)_3_1_Electric_Field_Lines"))
    end

    sim.detector = SolidStateDetector(sim.detector, ADLChargeDriftModel(T = T))

    calculate_drift_fields!(sim)

    pos = if key == :InvertedCoax
        CylindricalPoint{T}[ CylindricalPoint{T}( 0.02, deg2rad(10), 0.025 ) ]
    elseif key == :CGD
        CartesianPoint{T}[ CartesianPoint{T}( 0.006, 0.005, 0.005  ) ]
    elseif key == :BEGe
        CylindricalPoint{T}[ CylindricalPoint{T}( 0.016, deg2rad(10), 0.015  ) ]
    elseif key == :Coax
        CylindricalPoint{T}[ CylindricalPoint{T}( 0.016, deg2rad(10), 0.005  ) ]
    elseif key == :Spherical
        CylindricalPoint{T}[ CylindricalPoint{T}( 0.00, deg2rad(0), 0.0  ) ]
    end
    energy_depos = T[1460]
    @assert in(pos[1], sim.detector) "Test point $(pos[1]) not inside the detector $(key)."

    event = Event(pos, energy_depos)
    simulate!(event, sim)

    plot(sim.detector)
    plot!(event.drift_paths)
    savefig(joinpath(outputdir, "$(key)_4_charge_drift"))

    plot([event.waveforms...], size = (1200, 600), lw = 1.5)
    savefig(joinpath(outputdir, "$(key)_5_induced_signals"))

end

@info "Finished testing."
@info "Test output saved in: $outputdir"
