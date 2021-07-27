using BenchmarkTools
using SolidStateDetectors
using Plots

### Inner loop of SOR 
begin
    T = Float32
    sim = Simulation{T}(SSD_examples[:Coax]);

    imp_density = SolidStateDetectors.CylindricalImpurityDensity{Float32}((0.0f0, 0.0f0, 1f16), (0.0f0, 0.0f0, 0.0f0));
    sim.detector = SolidStateDetector(sim.detector, imp_density); 
    sim.detector = SolidStateDetector(sim.detector, contact_id = 1, contact_potential = 2000); # Make it undepleted

    calculate_electric_potential!( sim, depletion_handling = true,
                                        convergence_limit = 1e-4,
                                        refinement_limits = [0.2, 0.1, 0.05])

    plot(sim.point_types, φ = 30)

    fssrb = SolidStateDetectors.PotentialSimulationSetupRB( sim.detector, 
                                                            sim.electric_potential.grid, 
                                                            SolidStateDetectors.material_properties[SolidStateDetectors.materials["vacuum"]],
                                                            sim.electric_potential.data, 
                                                            sor_consts = (T(1), T(1)) );


    nthreads = 1 #Base.Threads.nthreads()
    is_weighting_potential = Val{false}()
    depletion_handling = Val{false}()
    only2d = Val{false}()
    update_even_points = Val{true}()
    update_uneven_points = Val{false}()
end

begin
    depletion_handling = Val{false}()
    @btime SolidStateDetectors.update!($fssrb, $nthreads, $update_even_points, $depletion_handling, $is_weighting_potential, $only2d)
end
begin
    depletion_handling = Val{true}()
    @btime SolidStateDetectors.update!($fssrb, $nthreads, $update_even_points, $depletion_handling, $is_weighting_potential, $only2d)
end


# Innerloop:
begin
    gw1 = fssrb.geom_weights[1].weights  # r or x 
    gw2 = fssrb.geom_weights[2].weights  # φ or y
    gw3 = fssrb.geom_weights[3].weights  # z or z
    
    rb_tar_idx = 1
    rb_src_idx = 2
    idx3 = 2

    depletion_handling = Val{false}()
    @btime SolidStateDetectors.innerloops!( $idx3, $rb_tar_idx, $rb_src_idx, 
            $gw1, $gw2, $gw3, $fssrb, $update_even_points, 
            $depletion_handling, $is_weighting_potential, $only2d)
    
    depletion_handling = Val{true}()
    @btime SolidStateDetectors.innerloops!( $idx3, $rb_tar_idx, $rb_src_idx, 
            $gw1, $gw2, $gw3, $fssrb, $update_even_points, 
            $depletion_handling, $is_weighting_potential, $only2d)
end

# apply_boundary_conditions!
begin 
    @btime SolidStateDetectors.apply_boundary_conditions!($fssrb, $update_even_points, $only2d)
    @btime SolidStateDetectors.apply_boundary_conditions!($fssrb, $update_uneven_points, $only2d)
end

begin
    depletion_handling = Val{false}()
    SolidStateDetectors.update!(fssrb; use_nthreads = nthreads, depletion_handling = depletion_handling,
        only2d = only2d, is_weighting_potential= is_weighting_potential)
    @btime SolidStateDetectors.update!($fssrb; use_nthreads = $nthreads, depletion_handling = $depletion_handling,
        only2d = $only2d, is_weighting_potential= $is_weighting_potential)
  
    depletion_handling = Val{true}()
    SolidStateDetectors.update!(fssrb; use_nthreads = nthreads, depletion_handling = depletion_handling,
        only2d = only2d, is_weighting_potential= is_weighting_potential)
    @btime SolidStateDetectors.update!($fssrb; use_nthreads = $nthreads, depletion_handling = $depletion_handling,
        only2d = $only2d, is_weighting_potential= $is_weighting_potential)
end