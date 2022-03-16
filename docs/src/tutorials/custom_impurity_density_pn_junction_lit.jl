# # Advanced Example: Custom Impurity Profile

# A very important input to the simulation is the impurity profile of the semiconductor.
# It influences the electric potential and, thus, the electric field, capacitances,
# drift paths and also the induced waveforms.

# Arbitrary impurity density profiles can be defined and assigned to the semiconductor of the detector.
# Also have a look into the manual: [Impurity Densities](@ref) and [Custom Impurity Density](@ref).

# As an example, we define the impurity profile of a simple [p-n junction](https://en.wikipedia.org/wiki/P%E2%80%93n_junction)
# and calculate the electric potential and depleted region and reproduce the 
# typical [plot](https://en.wikipedia.org/wiki/P%E2%80%93n_junction#/media/File:Pn-junction-equilibrium-graphs.png)
# shown in various books explaining p-n junctions.

# For this, we use one of the example detectors: the infinite parallel plate capacitor in cartesian coordinates.
# The plates span over `y` and `z`. It is made infinite by defining the world smaller than the contacts
# and using reflecting boundaries in `y` and `z`.
# Let's get started by loading the packages and the detector and have a look at its geometry:

using Plots
using SolidStateDetectors
using Unitful

sim = Simulation{Float32}(SSD_examples[:InfiniteParallelPlateCapacitor])
plot(sim.detector)

# ## Define the Impurity Density

# Now, we define the model for impurity density of the p-n junction.
# In order to do so, we need two define two things.

# 1. We have to define a `struct` which has to be 
# subtype of `SolidStateDetectors.AbstractImpurityDensity{T}`.
# We define the model with two parameters in order to be able to vary the density level 
# on the p-type and n-type regions independently. 
# More parameters could be added. E.g., the position of the p-n junction.

struct PNJunctionImpurityDensity{T} <: SolidStateDetectors.AbstractImpurityDensity{T} 
    p_type_density::T
    n_type_density::T
end

# 2. We have to define a method for `SolidStateDetectors.get_impurity_density` for our 
# just defined impurity density. 
# It has to take two arguments. An instance of our model and a point.
# Since the infinite parallel plate capacitor is configured to be calculated on a cartesian grid,
# the point has to be a `CartesianPoint`.
# The function returns the impurity density at the respective position.
# The returned value has to be in units of 1/m``^{3}``.

#=
!!! note "Sign of the impurity density"
    The sign of the impurity density determines whether it is p-type or n-type.

    p-type region <-> negative sign: Holes are the majority carriers and are free to move and diffuse into the n-type region. 
        Electrons are fixed in the lattice. Thus, a negative fixed space charge density is left behind in the depleted p-type region.

    n-type region <-> positive sign: Electrons are the majority carriers and are free to move and diffuse into the p-type region. 
        Holes are fixed in the lattice. Thus, a positive fixed space charge density is left behind in the depleted n-type region.
=#

function SolidStateDetectors.get_impurity_density(
    cdm::PNJunctionImpurityDensity{T}, 
    pt::CartesianPoint{T}
)::T where {T}
    x = pt[1] # In this example, we only need the `x` coordinate of the point.
    if x > 0
        -cdm.p_type_density # p-type region -> electrons are fixed -> negative charge 
    else 
        cdm.n_type_density  # n-type region -> holes are fixed -> positive charge
    end
end

# ## Assign the Density and Calculate the Fields

# Now, we create an instance of our model with reasonable parameters 
# (for the geometry and bias voltage of the example detector)
# and assign it to the detector.

pn_junction_impurity_density = PNJunctionImpurityDensity{Float32}(3e16, 1.5e16)

sim.detector = SolidStateDetector(sim.detector, pn_junction_impurity_density);

# We are now ready to calculate the electric potential. 
# Depletion handling has to be turned on, via the keyword `depletion_handling = true`,
# in order to detect the regions where the detector is depleted.
# The other settings for `calculate_electric_potential!` are very specific for this effectively 1D problem
# to be calculated in 3D.
# Afterwards, we also calculate the electric field and plot the 
# electric potential, impurity scale map and the point type map. 
# The two latter ons show us where the detector is depleted (`impurity scale == 1`).

calculate_electric_potential!(
    sim, 
    convergence_limit = 1e-6,
    max_tick_distance = (0.002, 1, 1) .* u"cm",
    min_tick_distance = (1e-7, 1, 1) .* u"m",
    refinement_limits = [0.005, 0.001],
    depletion_handling = true,
    verbose = true
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

# ## 1D Plots of the P-N Junction

# As this is in principle a 1D simulation, lets better plot the 
# charge density from the impurities, the electric potential and the electric field only over `x`.

xs = sim.electric_potential.grid[1] * 1000;
plot(
    begin
        plot(xs, sim.q_eff_imp.data[:, 1, 1] .* sim.imp_scale[:, 1, 1], lw = 4, label ="", yguide = "\$\\rho\$ [e/m\$^3\$]")
        vline!([sim.detector.contacts[1].geometry.origin.x] * 1000, lw = 4, color = "green", label = "N+ Contact")
        vline!([sim.detector.contacts[2].geometry.origin.x] * 1000, lw = 4, color = "red", label = "P+ Contact")
        vline!([0], lw = 4, color = "black", label = "p-n junction")
    end,
    begin
        plot(xs, sim.electric_potential.data[:,1,1], lw = 4, label ="", yguide = "\$\\Phi\$ [V]")
        vline!([sim.detector.contacts[1].geometry.origin.x] * 1000, lw = 4, color = "green", label = "N+ Contact")
        vline!([sim.detector.contacts[2].geometry.origin.x] * 1000, lw = 4, color = "red", label = "P+ Contact")
        vline!([0], lw = 4, color = "black", label = "p-n junction")
    end,
    begin
        plot(xs, map(E -> E[1]/1000, sim.electric_field.data[:,1,1]), lw = 4, label ="", yguide = "\$\\mathcal{E}_z\$ [V/mm]")
        vline!([sim.detector.contacts[1].geometry.origin.x] * 1000, lw = 4, color = "green", label = "N+ Contact")
        vline!([sim.detector.contacts[2].geometry.origin.x] * 1000, lw = 4, color = "red", label = "P+ Contact")
        vline!([0], lw = 4, color = "black", label = "p-n junction")
    end,
    layout = (3, 1),
    xlims = (-5.1, 5.1),
    size = (800, 800),
    lw = 2,
    xguide = "x [mm]"
)