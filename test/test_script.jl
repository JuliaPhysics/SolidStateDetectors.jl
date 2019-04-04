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
sleep(2)

plot() # creates a plot so that the plots during the following loop pop up.

for key in keys(SSD_examples)
    @info "Now test detector type: $key"
    det = SolidStateDetector{T}(SSD_examples[key])

    g = Grid(det)

    Init_setup = SolidStateDetectors.PotentialSimulationSetup(det);
    plot(Init_setup, size = (1200, 1200))
    savefig(joinpath(outputdir, "$(key)_init_setup"))
    
    for nrefs in [0, 1, 2]
        setup = calculate_electric_potential(det, max_refinements = nrefs)
        plot(setup, size = (1200, 1200))
        savefig(joinpath(outputdir, "$(key)_setup_$(nrefs)_refinements"))
    end    

    setup = calculate_electric_potential(det, max_refinements = 3 );
    plot(setup, size = (1200, 1200))
    savefig(joinpath(outputdir, "$(key)_setup_$(1)_refinements"))

    E_pot, point_types = if typeof(det) <: SolidStateDetector{<:AbstractFloat, :Cylindrical} && length(setup.grid[:φ]) == 1 # 2D 
        ElectricPotential(setup, n_points_in_φ = 36),
        PointTypes(setup,  n_points_in_φ = 36);
    else
        ElectricPotential(setup),
        PointTypes(setup);
    end;
    plot(E_pot, y = 0)

    E_field = SolidStateDetectors.get_electric_field_from_potential(E_pot, point_types);
    SSD.write_grid_to_detector!(det, E_pot.grid);

    # drift_model = ADLChargeDriftModel();
    drift_model = VacuumChargeDriftModel();
    electron_drift_field = get_electron_drift_field(E_field, drift_model);
    hole_drift_field = get_hole_drift_field(E_field, drift_model);

    electron_drift_field_interpolated = SolidStateDetectors.get_interpolated_drift_field(electron_drift_field, E_pot.grid);
    hole_drift_field_interpolated = SolidStateDetectors.get_interpolated_drift_field(hole_drift_field, E_pot.grid);

    pos = CartesianPoint{T}[ CartesianPoint{T}( 0.006, 0.000, 0.005  ) ] # this point should be inside all test detectors
    @assert in(pos[1], det) "Test point $(pos[1]) not inside the detector $(key)."

    drift_paths = SSD.drift_charges(det, pos, electron_drift_field_interpolated, hole_drift_field_interpolated);

end

@info "Finished testing."
@info "Test output saved in: $outputdir"


key = :CGD
det = SolidStateDetector{T}(SSD_examples[key])

setup = calculate_electric_potential(det, max_refinements = 3 );
plot(setup, size = (1200, 1200))
savefig(joinpath(outputdir, "$(key)_setup_$(1)_refinements"))

E_pot, point_types = if typeof(det) <: SolidStateDetector{<:AbstractFloat, :Cylindrical} && length(setup.grid[:φ]) == 1 # 2D 
    ElectricPotential(setup, n_points_in_φ = 36),
    PointTypes(setup,  n_points_in_φ = 36);
else
    ElectricPotential(setup),
    PointTypes(setup);
end;
plot(E_pot, y = 0)

plot(E_pot.grid[:x], E_pot.data[:, 120, 120])


E_field = SolidStateDetectors.get_electric_field_from_potential(E_pot, point_types);
SSD.write_grid_to_detector!(det, E_pot.grid);

# drift_model = ADLChargeDriftModel();
drift_model = VacuumChargeDriftModel();
electron_drift_field = get_electron_drift_field(E_field, drift_model);
hole_drift_field = get_hole_drift_field(E_field, drift_model);

electron_drift_field_interpolated = SolidStateDetectors.get_interpolated_drift_field(electron_drift_field, E_pot.grid);
hole_drift_field_interpolated = SolidStateDetectors.get_interpolated_drift_field(hole_drift_field, E_pot.grid);

pos = CartesianPoint{T}[ CartesianPoint{T}( 0,0,0)]#0.006, 0.000, 0.005  ) ] # this point should be inside all test detectors
@assert in(pos[1], det) "Test point $(pos[1]) not inside the detector $(key)."

drift_paths = SSD.drift_charges(det, pos, electron_drift_field_interpolated, hole_drift_field_interpolated);

plot(det.crystal_geometry.a, lc=:black, lw = 2, label="N-type bulk")

plot!(drift_paths)

plot!(drift_paths[1].e_path[1:2], lw= 2)
plot!(drift_paths[1].e_path[1:220], lw= 2)


plot!(det.contacts[1].geometry[1], lc=:orange, lw = 2, label="N+")
plot!(det.contacts[2].geometry[1], lc=:blue, lw = 2, label="P+")
plot!(drift_paths, xlabel = "x", ylabel = "y", zlabel = "z")

efs = SSD.get_electric_field_strength(E_field) 

plot(efs[:, 100, 100])

plot(efs[:, 95, :]', st=:heatmap, clims = (0, 2000 / 0.02 * 2))