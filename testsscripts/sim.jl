using Plots
using SolidStateDetectors
using Unitful
import LegendHDF5IO

T = Float32

sim_file = "$(@__DIR__)/ssd-sims.h5"
sim = if ispath(sim_file)
    ssd_read(sim_file, Simulation)
else
    sim = Simulation{T}(SSD_examples[:InvertedCoax])

    plot(sim.detector, xunit = u"mm", yunit = u"mm", zunit = u"mm", size = (700, 700))

    calculate_electric_potential!( sim, refinement_limits = [0.2, 0.1, 0.05, 0.01])

    calculate_electric_field!(sim, n_points_in_φ = 72)


    charge_drift_model = ADLChargeDriftModel()
    sim.detector = SolidStateDetector(sim.detector, charge_drift_model)

    ssd_write(sim_file, sim)
    sim
end

#=
plot(
    plot(sim.electric_potential, φ = 20), # initial electric potential (boundary conditions)
    plot(sim.point_types), # map of different point types: fixed point / inside or outside detector volume / depleted/undepleted
    plot(sim.q_eff_imp), # charge density distribution
    plot(sim.ϵ_r), # dielectric distribution
    layout = (1, 4), size = (1600, 500)
)

plot(
    sim.electric_field, full_det = true, φ = 0.0, size = (700, 700),
    xunit = u"mm", yunit = u"mm", zunit = u"V/mm", clims = (0,800).*u"V/mm"
)

plot_electric_fieldlines!(sim, full_det = true, φ = 0.0)
=#