using Plots
using SolidStateDetectors
using Unitful

sim = Simulation{Float32}(SSD_examples[:InfiniteParallelPlateCapacitor])
plot( sim.detector.semiconductor, fillalpha = 0.4)
plot!(sim.detector, xlims = (-0.006, 0.006), aspect_ratio = :none)

struct PNJunctionImpurityDensity{T} <: SolidStateDetectors.AbstractImpurityDensity{T}
    p_type_density::T
    n_type_density::T
end

function SolidStateDetectors.get_impurity_density(
    cdm::PNJunctionImpurityDensity{T},
    pt::SolidStateDetectors.AbstractCoordinatePoint{T}
)::T where {T}
    cpt = CartesianPoint(pt)
    x = cpt[1] # In this example, we only need the `x` coordinate of the point.
    if x > 0
        -cdm.p_type_density # p-type region -> electrons are fixed -> negative charge
    else
        cdm.n_type_density  # n-type region -> holes are fixed -> positive charge
    end
end

pn_junction_impurity_density = PNJunctionImpurityDensity{Float32}(3e16, 1.5e16)
sim.detector = SolidStateDetector(sim.detector, pn_junction_impurity_density);

calculate_electric_potential!(
    sim,
    convergence_limit = 1e-6,
    max_tick_distance = (0.002, 1, 1) .* u"cm",
    min_tick_distance = (1e-7, 1, 1) .* u"m",
    refinement_limits = [0.005, 0.001],
    depletion_handling = true
)
calculate_electric_field!(sim)

plot(
    plot(sim.electric_potential, y = 0),
    plot(sim.imp_scale, y = 0),
    plot(sim.point_types, y = 0),
    layout = (3, 1),
    aspect_ratio = :none,
    size = (800, 800)
)

xs = uconvert.(u"mm", sim.electric_potential.grid[1] * u"m");
plot(
    begin
        ρ_x = map(x -> SolidStateDetectors.get_impurity_density(
            sim.detector.semiconductor.impurity_density_model, CartesianPoint{Float32}(x, 0, 0)),
            sim.electric_potential.grid[1]
        )
        plot(xs, ρ_x .* sim.imp_scale[:, 1, 1], lw = 4, label ="", xguide = "x", unitformat = :square, yguide = "\$\\rho\$ [e/m\$^3\$]")
        vline!([sim.detector.contacts[1].geometry.origin.x] * 1000, lw = 4, color = "green", label = "N+ Contact")
        vline!([sim.detector.contacts[2].geometry.origin.x] * 1000, lw = 4, color = "red", label = "P+ Contact")
        vline!([0], lw = 4, color = "black", label = "p-n junction")
    end,
    begin
        plot(xs, sim.electric_potential.data[:,1,1], lw = 4, label ="", xguide = "x", unitformat = :square, yguide = "\$\\Phi\$ [V]")
        vline!([sim.detector.contacts[1].geometry.origin.x] * 1000, lw = 4, color = "green", label = "N+ Contact")
        vline!([sim.detector.contacts[2].geometry.origin.x] * 1000, lw = 4, color = "red", label = "P+ Contact")
        vline!([0], lw = 4, color = "black", label = "p-n junction")
    end,
    begin
        plot(xs, map(E -> E[1]/1000, sim.electric_field.data[:,1,1]), lw = 4, label ="", xguide = "x", unitformat = :square, yguide = "\$\\mathcal{E}_x\$ [V/mm]")
        vline!([sim.detector.contacts[1].geometry.origin.x] * 1000, lw = 4, color = "green", label = "N+ Contact")
        vline!([sim.detector.contacts[2].geometry.origin.x] * 1000, lw = 4, color = "red", label = "P+ Contact")
        vline!([0], lw = 4, color = "black", label = "p-n junction")
    end,
    layout = (3, 1),
    xlims = (-5.1, 5.1),
    size = (800, 800),
    lw = 2,
)
