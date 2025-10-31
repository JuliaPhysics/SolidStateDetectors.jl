var documenterSearchIndex = {"docs": [

{
    "location": "#",
    "page": "Home",
    "title": "Home",
    "category": "page",
    "text": ""
},

{
    "location": "#SolidStateDetectors.jl-1",
    "page": "Home",
    "title": "SolidStateDetectors.jl",
    "category": "section",
    "text": "SolidStateDetectors.jl is a Julia package for fast 2D and 3D simulation of Solid State Detectors.Pages = [\"man/installation.md\"]Pages = [\"man/detectors.md\"]Pages = [\"man/electric_potentials.md\"]Pages = [\"man/weighting_potentials.md\"]Pages = [\"man/electric_fields.md\"]Pages = [\"man/drift_fields.md\"]Pages = [\"man/IO.md\"]"
},

{
    "location": "man/installation/#",
    "page": "Installation",
    "title": "Installation",
    "category": "page",
    "text": ""
},

{
    "location": "man/installation/#Installation-1",
    "page": "Installation",
    "title": "Installation",
    "category": "section",
    "text": "This package is not a registered package (yet).Install viausing Pkg; pkg\"add https://github.com/JuliaHEP/SolidStateDetectors.jl.git\""
},

{
    "location": "man/installation/#Vizualization-/-Plotting-(Optional)-1",
    "page": "Installation",
    "title": "Vizualization / Plotting (Optional)",
    "category": "section",
    "text": "This package provides serveral plot recipes for different outputs for the plotting package Plots.jl.In order to use these also install the Plots.jl package viausing Pkg; pkg\"add Plots\"It is recommended to use pyplot as backend. Install viausing Pkg; pkg\"add PyPlot\"Then you can load it viausing Plots; pyplot();"
},

{
    "location": "man/detectors/#",
    "page": "Detectors",
    "title": "Detectors",
    "category": "page",
    "text": ""
},

{
    "location": "man/detectors/#Detectors-1",
    "page": "Detectors",
    "title": "Detectors",
    "category": "section",
    "text": "Currently, three classes of detectors are supported: Coaxial, Inverted Coax, and BEGe type detectors. All detector properties are specified in a one-for-all .json file. "
},

{
    "location": "man/detectors/#Example-1)-Inverted-Coax-1",
    "page": "Detectors",
    "title": "Example 1) Inverted Coax",
    "category": "section",
    "text": "Example minimum config file for an Inverted Coax detector (IVC) plus explanations.  Remember, comments are not allowed in JSON files and have to be deleted if you want to use it.{\n    \"name\":\"ExampleInvertedCoax\", // Arbitrary name of the detector\n    \"class\":\"InvertedCoax\",  // either \"Coax\", \"BEGe\", \"InvertedCoax\"\n    \"type\":\"p\", // either \"p\", \"ptype\", \"p-type\", \"n\", \"ntype\" or \"n-type\"\n    \"cyclic\":0, // The periodicity of the detector in degree. \n                // `0` means complete symmetric in θ -> 2D simulation. \n                // The condition \"360 / `cyclic` = n; n = 1, 2, 3, ...\" must be fulfilled.\n    \"mirror_symmetry_θ\": \"true\", // set to true, if a mirror symmetry exists within the periodicity specified via `cyclic`\n    \"materials\":{\n        \"detector\":\"HPGe\",  // Material of the detector \n        \"environment\":\"Vacuum\" // Material of the environment\n    },\n    \"geometry\":\n    {\n        \"unit\": \"mm\", // unit of the values which the user specifies in this config file\n        \"crystal\":\n        {\n            \"length\": 80.0, // crystal goes from z=0mm to z=80mm\n            \"radius\": 35.0  // crystal goes from r=0mm to r=35mm\n        },\n        \"point_contact\":\n        {\n            \"endplate\":\"bot\",\n            \"depth\":3e-4,\n            \"radius\":9.5\n        },\n        \"groove\":\n        {\n            \"endplate\":\"bot\",\n            \"rInner\":13.5,\n            \"width\":0.0,\n            \"depth\":0.0\n        },\n        \"borehole\": // IVC: Borehole starts at the top of the crystal \n        {\n          \"length\":55.0,\n          \"radius\":5.0\n        },\n        \"taper\": \n        {\n            \"outer\":\n            {\n                \"rInner\":24.42,\n                \"length\":60\n            },\n            \"inner\":\n            {\n                \"rOuter\":0,\n                \"length\":0\n            }\n        }\n    },\n    \"charge_carrier_density\":\n    {\n        \"unit\": \"cm-3\",\n        \"top\": 0.94, // times 10^10\n        \"bot\": 1.10  // just the absolute value. The sign is determined through the bulk material\n    },\n    \"segmentation\":\n    {\n        \"n_contacts_total\": 3, // n-contacts + p-contacts\n        \"core\": // core = n-contact of the detector\n        {\n            \"type\": \"Tubs\", // \"Box\" volume in cylindrical coordinates\n            \"rStart\": 0.0,\n            \"rStop\": 3.0,\n            \"phiStart\": 0,\n            \"phiStop\": 360,\n            \"zStart\": 0.0,\n            \"zStop\": 0.0,\n            \"potential\": 0.0 // Bias voltage of the core\n        },        \n        \"n_individual_segments\": 2, // now the p-contacts\n        \"S1\":\n        {\n            \"type\": \"Tubs\",\n            \"rStart\": 15.0,\n            \"rStop\": 35.0,\n            \"phiStart\": 0,\n            \"phiStop\": 360,\n            \"zStart\": 0,\n            \"zStop\": 0,\n            \"potential\": 3500.0, // Bias voltage of this segment\n            \"boundaryWidth\":\n            {\n                \"radial\":0.0,\n                \"vertical\": 0.0,\n                \"horizontal\": 0.0\n            },\n            \"repetitive\":false,\n            \"repetitions\":\n            {\n                \"radial\":0,\n                \"vertical\": 0,\n                \"horizontal\": 5\n            }\n        },\n        \"S2\":\n        {\n            \"type\": \"Tubs\",\n            \"rStart\": 35.0,\n            \"rStop\": 35.0,\n            \"phiStart\": 0,\n            \"phiStop\": 360,\n            \"zStart\": 0.0,\n            \"zStop\": 20.0,\n            \"potential\": 3500.0, // Bias voltage of this segment\n            \"boundaryWidth\":\n            {\n                \"radial\":0.0,\n                \"vertical\": 0.0,\n                \"horizontal\": 0.0\n            },\n            \"repetitive\":false,\n            \"repetitions\":\n            {\n                \"radial\": 0,\n                \"vertical\": 0,\n                \"horizontal\": 5\n            }\n        },\n        \"custom_grouping\":false, // Multiple geometrical segments can build up one contact:\n                                 // Specify which geometrical segments belong to which contact\n                                 // To specify geometrical segments `3`, `4`, `5`, `6` use \"3-4-5-6\" (not only \"3-6\")\n        \"Chn1\":\"1\",    // contact \"Chn1\" is build through the geometrical segment `1`. (This is the core in this case)\n        \"Chn2\":\"2-3\"   // contact \"Chn2\" is build through the geometrical segments `2` & `3`. \n    }\n}"
},

{
    "location": "man/electric_potentials/#",
    "page": "Electric Potentials",
    "title": "Electric Potentials",
    "category": "page",
    "text": ""
},

{
    "location": "man/electric_potentials/#Electric-Potentials-1",
    "page": "Electric Potentials",
    "title": "Electric Potentials",
    "category": "section",
    "text": ""
},

{
    "location": "man/electric_potentials/#Simulation-Algorithm-1",
    "page": "Electric Potentials",
    "title": "Simulation Algorithm",
    "category": "section",
    "text": "The electric potential is calculated through successive over relaxation (SOR).The calculation is based on Gauss\' law in matternabla epsilon_r(vecr) nabla varphi(vecr) = dfracrho(vecr)epsilon_0which is numerically solved on a 3-dimensional adaptive red/black grid.The red/black division allows for multithreading and the adaptive grid saves computation time since it only increases the grid point density in areas where it is critical.For now, only cylindrical coordinates are supported. "
},

{
    "location": "man/electric_potentials/#Multithreading-1",
    "page": "Electric Potentials",
    "title": "Multithreading",
    "category": "section",
    "text": "To use multiple threads for the simulation, the environement variable JULIA_NUM_THREADS must be set before Julia is started. In case of bash this is done throughexport JULIA_NUM_THREADS=4Note that the user still has to parse the number of threads to the function as a keyword nthreads. See calculate_electric_potential and calculate_weighting_potential."
},

{
    "location": "man/electric_potentials/#Example-1",
    "page": "Electric Potentials",
    "title": "Example",
    "category": "section",
    "text": "To simulate a detector the detector must first be defined (SolidStateDetector).using Plots; pyplot()\nusing SolidStateDetectors\n\ndetector = SolidStateDetector(SSD_examples[:InvertedCoax]);Then, the electric potential can be calculated through the function calculate_electric_potential. It returns a collection struct SolidStateDetectors.PotentialSimulationSetup in which the electric potential is stored. Also it stores the grid, the charge density, the dielectric distribution and the point types PointTypes. The electric potential can be extracted via the function ElectricPotential(::SolidStateDetectors.PotentialSimulationSetup). The struct ElectricPotential also holds the corresponding grid Grid. E_pot_setup = calculate_electric_potential(detector);\nE_pot = ElectricPotential(E_pot_setup, n_points_in_θ = 36);A plot recipe for the struct ElectricPotential exists so the result can be visualized throughplot(E_pot, θ=0)It can also be plotted in the r-θ plane:plot(E_pot, z=0.04)"
},

{
    "location": "man/weighting_potentials/#",
    "page": "Weighting Potentials",
    "title": "Weighting Potentials",
    "category": "page",
    "text": ""
},

{
    "location": "man/weighting_potentials/#Weighting-Potentials-1",
    "page": "Weighting Potentials",
    "title": "Weighting Potentials",
    "category": "section",
    "text": ""
},

{
    "location": "man/weighting_potentials/#Simulation-Algorithm-1",
    "page": "Weighting Potentials",
    "title": "Simulation Algorithm",
    "category": "section",
    "text": "The weighting potential for an electrode is internally calculated with the same function as the electric potential.The only difference is that the charge carrier density ρ(vecr) is set to 0 and all electrodes are fixed to 0 but the one electrode which is fixed to 1."
},

{
    "location": "man/electric_fields/#",
    "page": "Electric Fields",
    "title": "Electric Fields",
    "category": "page",
    "text": ""
},

{
    "location": "man/electric_fields/#Electric-Fields-1",
    "page": "Electric Fields",
    "title": "Electric Fields",
    "category": "section",
    "text": ""
},

{
    "location": "man/drift_fields/#",
    "page": "Drift Fields",
    "title": "Drift Fields",
    "category": "page",
    "text": ""
},

{
    "location": "man/drift_fields/#Drift-Fields-1",
    "page": "Drift Fields",
    "title": "Drift Fields",
    "category": "section",
    "text": ""
},

{
    "location": "man/IO/#",
    "page": "IO",
    "title": "IO",
    "category": "page",
    "text": ""
},

{
    "location": "man/IO/#IO-1",
    "page": "IO",
    "title": "IO",
    "category": "section",
    "text": ""
},

{
    "location": "api/#",
    "page": "API",
    "title": "API",
    "category": "page",
    "text": ""
},

{
    "location": "api/#API-1",
    "page": "API",
    "title": "API",
    "category": "section",
    "text": "DocTestSetup  = quote\n    using SolidStateDetectors\nend"
},

{
    "location": "api/#Types-1",
    "page": "API",
    "title": "Types",
    "category": "section",
    "text": "Order = [:type]"
},

{
    "location": "api/#Functions-1",
    "page": "API",
    "title": "Functions",
    "category": "section",
    "text": "Order = [:function]"
},

{
    "location": "api/#SolidStateDetectors.ADLChargeDriftModel",
    "page": "API",
    "title": "SolidStateDetectors.ADLChargeDriftModel",
    "category": "type",
    "text": "ADLChargeDriftModel{T <: AbstractFloat} <: AbstractChargeDriftModels\n\nFields\n\nelectrons::CarrierParameters{T}\nholes::CarrierParameters{T}\n`phi110::T\ngammas::SVector{4, SMatrix{3,3,T}}\n\n\n\n\n\n"
},

{
    "location": "api/#SolidStateDetectors.ElectricPotential-Union{Tuple{PotentialSimulationSetup{T,3,:Cylindrical}}, Tuple{T}} where T",
    "page": "API",
    "title": "SolidStateDetectors.ElectricPotential",
    "category": "method",
    "text": "ElectricPotential(setup::PotentialSimulationSetup{T, 3, :Cylindrical} ; kwargs...)::ElectricPotential{T, 3, :Cylindrical}\n\nExtracts the electric potential from setup and extrapolate it to an 2π grid.\n\nFor 2D grids (r and z) the user has to set the keyword n_points_in_θ::Int, e.g.: n_points_in_θ = 36.\n\n\n\n\n\n"
},

{
    "location": "api/#SolidStateDetectors.ElectricPotential-Union{Tuple{SolidStateDetector{T}}, Tuple{T}} where T",
    "page": "API",
    "title": "SolidStateDetectors.ElectricPotential",
    "category": "method",
    "text": "ElectricPotential(detector::SolidStateDetector{T}; kwargs...) where {T}\n\nCalls calculate_electric_potential(detector; kwargs...) end extracts only the electric potential and returns it.\n\n\n\n\n\n"
},

{
    "location": "api/#SolidStateDetectors.Grid",
    "page": "API",
    "title": "SolidStateDetectors.Grid",
    "category": "type",
    "text": "T: tick type\nN: N dimensional\nS: System (Cartesian, Cylindrical...)\n\n\n\n\n\n"
},

{
    "location": "api/#SolidStateDetectors.PointTypes",
    "page": "API",
    "title": "SolidStateDetectors.PointTypes",
    "category": "type",
    "text": "PointTypes{T, N, S} <: AbstractArray{T, N}\n\nPointTypes stores the point type of each grid point.\n\n\n\n\n\n"
},

{
    "location": "api/#SolidStateDetectors.SolidStateDetector",
    "page": "API",
    "title": "SolidStateDetectors.SolidStateDetector",
    "category": "type",
    "text": "abstract type SolidStateDetector{T}\n\nSupertype of all detector structs.\n\n\n\n\n\n"
},

{
    "location": "api/#SolidStateDetectors.SolidStateDetector-Union{Tuple{AbstractString}, Tuple{T}} where T<:AbstractFloat",
    "page": "API",
    "title": "SolidStateDetectors.SolidStateDetector",
    "category": "method",
    "text": "SolidStateDetector{T}(filename::AbstractString)::SolidStateDetector where T <: AbstractFloat\n\nReads in a config-JSON file and returns an Detector struct which holds all information specified in the config file.\n\n\n\n\n\n"
},

{
    "location": "api/#SolidStateDetectors.VacuumChargeDriftModel",
    "page": "API",
    "title": "SolidStateDetectors.VacuumChargeDriftModel",
    "category": "type",
    "text": "VacuumChargeDriftModel <: AbstractChargeDriftModels\n\n\n\n\n\n"
},

{
    "location": "api/#SolidStateDetectors.WeightingPotential-Union{Tuple{T}, Tuple{SolidStateDetector{T},Int64}} where T",
    "page": "API",
    "title": "SolidStateDetectors.WeightingPotential",
    "category": "method",
    "text": "WeightingPotential(detector::SolidStateDetector{T}, channel::Int; kwargs...) where {T}\n\nCalls calculate_weighting_potential(detector, channel; kwargs...) end extracts only the weighting potential and returns it.\n\n\n\n\n\n"
},

{
    "location": "api/#SolidStateDetectors.add_fano_noise",
    "page": "API",
    "title": "SolidStateDetectors.add_fano_noise",
    "category": "function",
    "text": "add_fano_noise(E_dep::RealQuantity, E_ionisation::RealQuantity, f_fano::Real)::RealQuantity\n\nAdd Fano noise to an energy deposition E_dep, assuming a detector material ionisation energy E_ionisation and a Fano factor f_fano.\n\n\n\n\n\n"
},

{
    "location": "api/#SolidStateDetectors.calculate_electric_potential-Union{Tuple{SolidStateDetector{T}}, Tuple{S}, Tuple{N}, Tuple{T}} where S where N where T",
    "page": "API",
    "title": "SolidStateDetectors.calculate_electric_potential",
    "category": "method",
    "text": "calculate_electric_potential(det::SolidStateDetector{T}; <keyword arguments>) where {T}\n\nCompute the electric potential for the given Detector det on an adaptive grid through successive over relaxation. It returns a collection struct PotentialSimulationSetup{T} which stores the potential, the charge density, the dielectric distribution, pointtypes and the final grid.\n\nThere are serveral <keyword arguments> which can be used to tune the computation:\n\nKeywords\n\nconvergence_limit::Real: convergence_limit times the bias voltage sets the convergence limit of the relaxation. The convergence value is the absolute maximum difference of the potential between two iterations of all grid points. Default of convergence_limit is 5e-6 (times bias voltage).\nmax_refinements::Int: number of maximum refinements. Default is 2. Set it to 0 to switch off refinement.\nrefinement_limits::Vector{Real}: vector of refinement limits for each dimension (in case of cylindrical coordinates the order is r, θ, z). A refinement limit (e.g. refinement_limits[1]) times the bias voltage of the detector det is the maximum allowed voltage difference between two neighbouring grid points in the respective dimension. When the difference is larger, new points are created inbetween. Default is [1e-4, 1e-4, 1e-4].\ninit_grid_spacing::Vector{Real}: vector of the initial distances between two grid points for each dimension. For normal coordinates the unit is meter. For angular coordinates, the unit is radiance. It prevents the refinement to make the grid to fine. Default is [0.005, 10.0, 0.005]`.\nmin_grid_spacing::Vector{Real}: vector of the mimimum allowed distance between two grid points for each dimension. For normal coordinates the unit is meter. For angular coordinates, the unit is radiance. It prevents the refinement to make the grid to fine. Default is [1e-4, 1e-2, 1e-4].\ngrid::Grid{T, N, S}: Initial grid used to start the simulation. Default is Grid(detector, init_grid_spacing=init_grid_spacing).\ndepletion_handling::Bool: enables the handling of undepleted regions. Default is false.\nuse_nthreads::Int: Number of threads to use in the computation. Default is Base.Threads.nthreads(). The environment variable JULIA_NUM_THREADS must be set appropriately before the Julia session was started (e.g. export JULIA_NUM_THREADS=8 in case of bash).\nsor_consts::Vector{<:Real}: Two element array. First element contains the SOR constant for r = 0. Second contains the constant at the outer most grid point in r. A linear scaling is applied in between. First element should be smaller than the second one and both should be ∈ [1.0, 2.0]. Default is [1.4, 1.85].\nmax_n_iterations::Int: Set the maximum number of iterations which are performed after each grid refinement. Default is 10000. If set to -1 there will be no limit.\nverbose::Bool=true: Boolean whether info output is produced or not.\ninit_grid_spacing::Vector{<:Real}: Initial spacing of the grid. Default is [2e-3, 5, 2e-3] <=> [2mm, 5 degree, 2mm ]\n\n\n\n\n\n"
},

{
    "location": "api/#SolidStateDetectors.calculate_weighting_potential-Union{Tuple{S}, Tuple{N}, Tuple{T}, Tuple{SolidStateDetector{T},Int64}} where S where N where T",
    "page": "API",
    "title": "SolidStateDetectors.calculate_weighting_potential",
    "category": "method",
    "text": "calculate_weighting_potential(det::SolidStateDetector{T}; <keyword arguments>) where {T}\n\nCompute the weighting potential for the given Detector det on an adaptive grid through successive over relaxation. It returns a collection struct PotentialSimulationSetup{T} which stores the potential, the charge density, the dielectric distribution, pointtypes and the final grid.\n\nThere are serveral <keyword arguments> which can be used to tune the computation:\n\nKeywords\n\nconvergence_limit::Real: convergence_limit times the bias voltage sets the convergence limit of the relaxation. The convergence value is the absolute maximum difference of the potential between two iterations of all grid points. Default of convergence_limit is 5e-6 (times bias voltage).\nmax_refinements::Int: number of maximum refinements. Default is 2. Set it to 0 to switch off refinement.\nrefinement_limits::Vector{Real}: vector of refinement limits for each dimension (in case of cylindrical coordinates the order is r, θ, z). A refinement limit (e.g. refinement_limits[1]) times the bias voltage of the detector det is the maximum allowed voltage difference between two neighbouring grid points in the respective dimension. When the difference is larger, new points are created inbetween. Default is [1e-4, 1e-4, 1e-4].\ninit_grid_spacing::Vector{Real}: vector of the initial distances between two grid points for each dimension. For normal coordinates the unit is meter. For angular coordinates, the unit is radiance. It prevents the refinement to make the grid to fine. Default is [0.005, 10.0, 0.005]`.\nmin_grid_spacing::Vector{Real}: vector of the mimimum allowed distance between two grid points for each dimension. For normal coordinates the unit is meter. For angular coordinates, the unit is radiance. It prevents the refinement to make the grid to fine. Default is [1e-4, 1e-2, 1e-4].\ngrid::Grid{T, N, S}: Initial grid used to start the simulation. Default is Grid(detector, init_grid_spacing=init_grid_spacing).\ndepletion_handling::Bool: enables the handling of undepleted regions. Default is false.\nuse_nthreads::Int: Number of threads to use in the computation. Default is Base.Threads.nthreads(). The environment variable JULIA_NUM_THREADS must be set appropriately before the Julia session was started (e.g. export JULIA_NUM_THREADS=8 in case of bash).\nsor_consts::Vector{<:Real}: Two element array. First element contains the SOR constant for r = 0. Second contains the constant at the outer most grid point in r. A linear scaling is applied in between. First element should be smaller than the second one and both should be ∈ [1.0, 2.0]. Default is [1.4, 1.85].\nmax_n_iterations::Int: Set the maximum number of iterations which are performed after each grid refinement. Default is 10000. If set to -1 there will be no limit.\nverbose::Bool=true: Boolean whether info output is produced or not.\ninit_grid_spacing::Vector{<:Real}: Initial spacing of the grid. Default is [2e-3, 5, 2e-3] <=> [2mm, 5 degree, 2mm ]\n\n\n\n\n\n"
},

{
    "location": "api/#SolidStateDetectors.get_active_volume-Union{Tuple{PointTypes{T,3,:Cylindrical}}, Tuple{T}} where T",
    "page": "API",
    "title": "SolidStateDetectors.get_active_volume",
    "category": "method",
    "text": "get_active_volume(pts::PointTypes{T}) where {T}\n\nReturns an approximation of the active volume of the detector by summing up the cell volumes of all depleted cells.\n\n\n\n\n\n"
},

{
    "location": "api/#SolidStateDetectors.AbstractConfig",
    "page": "API",
    "title": "SolidStateDetectors.AbstractConfig",
    "category": "type",
    "text": "abstract type AbstractConfig{T <: AbstractFloat} end\n\nT: Type of points or precision. \n\n\n\n\n\n"
},

{
    "location": "api/#SolidStateDetectors.DiscreteAxis",
    "page": "API",
    "title": "SolidStateDetectors.DiscreteAxis",
    "category": "type",
    "text": "DiscreteAxis{T, BL, BR} <: AbstractAxis{T, BL, BR}\n\nT: Type of ticks\nBL, BR ∈ {:periodic, :reflecting, :infinite, :r0} \nBL: left boundary condition\nBR: right boundary condition\n\n\n\n\n\n"
},

{
    "location": "api/#SolidStateDetectors.DiscreteAxis-Union{Tuple{T}, Tuple{T,T,Symbol,Symbol,Symbol,Symbol,AbstractArray{T,1}}} where T",
    "page": "API",
    "title": "SolidStateDetectors.DiscreteAxis",
    "category": "method",
    "text": "DiscreteAxis(left_endpoint::T, right_endpoint::T, BL::Symbol, BR::Symbol, L::Symbol, R::Symbol, ticks::AbstractVector{T}) where {T}\n\nT: Type of ticks\nBL, BR ∈ {:periodic, :reflecting, :infinite, :r0} \nL, R {:closed, :open} \nticks: Ticks of the axis\n\n\n\n\n\n"
},

{
    "location": "api/#SolidStateDetectors.PointType",
    "page": "API",
    "title": "SolidStateDetectors.PointType",
    "category": "type",
    "text": "const PointType = UInt8\n\nStores certain information about a grid point via bit-flags. \n\nRight now there are:\n\n`const update_bit      = 0x01`\n`const undepleted_bit  = 0x02`\n`const pn_junction_bit = 0x04`\n\nHow to get information out of a PointType variable pt:\n\npt & update_bit == 0 -> do not update this point (for fixed points)     \npt & update_bit >  0 -> do update this point    \npt & undepleted_bit > 0 -> this point is undepleted\npt & pn_junction_bit > 0 -> this point belongs to the solid state detector. So it is in the volume of the n-type or p-type material.\n\n\n\n\n\n"
},

{
    "location": "api/#SolidStateDetectors.PotentialSimulationSetup",
    "page": "API",
    "title": "SolidStateDetectors.PotentialSimulationSetup",
    "category": "type",
    "text": "PotentialSimulationSetup{T, N, S} <: AbstractPotentialSimulationSetup{T, N}\n\nCollection struct. It holds the grid, the potential, the point types, the charge density and the dielectric distribution.\n\n\n\n\n\n"
},

{
    "location": "api/#SolidStateDetectors.RBArray-Union{Tuple{N}, Tuple{T}, Tuple{Type,Grid{T,N,:Cylindrical}}} where N where T",
    "page": "API",
    "title": "SolidStateDetectors.RBArray",
    "category": "method",
    "text": "RBExtBy2Array( et::Type, g::Grid{T, N, :Cylindrical} )::Array{et, N + 1} where {T, N}\n\nReturns a RedBlack array for the grid g. \n\n\n\n\n\n"
},

{
    "location": "api/#SolidStateDetectors.RBExtBy2Array-Union{Tuple{N}, Tuple{T}, Tuple{Type,Grid{T,N,:Cylindrical}}} where N where T",
    "page": "API",
    "title": "SolidStateDetectors.RBExtBy2Array",
    "category": "method",
    "text": "RBExtBy2Array( et::Type, g::Grid{T, N, :Cylindrical} )::Array{et, N + 1} where {T, N}\n\nReturns a RedBlack array for the grid g. The RedBlack array is extended in its size by 2 in each geometrical dimension.\n\n\n\n\n\n"
},

{
    "location": "api/#SolidStateDetectors.get_path_to_example_config_files-Tuple{}",
    "page": "API",
    "title": "SolidStateDetectors.get_path_to_example_config_files",
    "category": "method",
    "text": "get_path_to_example_config_files()::String\n\nReturns the path to example config files provided by the package.\n\n\n\n\n\n"
},

{
    "location": "api/#SolidStateDetectors.get_rbidx_right_neighbour-Tuple{Int64,Val{true},Val{true}}",
    "page": "API",
    "title": "SolidStateDetectors.get_rbidx_right_neighbour",
    "category": "method",
    "text": "get_rbidx_right_neighbour(rbidx::Int, ::Val{true}, ::Val{true})::Int\n\nneeds docu...\n\n\n\n\n\n"
},

{
    "location": "api/#SolidStateDetectors.innerloops!-Union{Tuple{_is_weighting_potential}, Tuple{_bulk_is_ptype}, Tuple{depletion_handling_enabled}, Tuple{even_points}, Tuple{T}, Tuple{Int64,Int64,Int64,Array{T,2},Array{T,2},Array{T,2},PotentialSimulationSetupRB{T,3,4,:Cylindrical},Val{even_points},Val{depletion_handling_enabled},Val{_bulk_is_ptype},Val{_is_weighting_potential}}} where _is_weighting_potential where _bulk_is_ptype where depletion_handling_enabled where even_points where T",
    "page": "API",
    "title": "SolidStateDetectors.innerloops!",
    "category": "method",
    "text": "innerloops!(  ir::Int, rb_tar_idx::Int, rb_src_idx::Int, gw_r::Array{T, 2}, gw_θ::Array{T, 2}, gw_z::Array{T, 2}, fssrb::PotentialSimulationSetupRB{T, 3, 4, :Cylindrical},\n                            update_even_points::Val{even_points},\n                            depletion_handling::Val{depletion_handling_enabled},\n                            bulk_is_ptype::Val{_bulk_is_ptype}  )::Nothing where {T, even_points, depletion_handling_enabled, _bulk_is_ptype}\n\n(Vectorized) inner loop for Cylindrical coordinates. This function does all the work in the fied calculation.                            \n\n\n\n\n\n"
},

{
    "location": "api/#SolidStateDetectors.nidx-Tuple{Int64,Val{true},Val{true}}",
    "page": "API",
    "title": "SolidStateDetectors.nidx",
    "category": "method",
    "text": "nidx( rbidx::Int, ::Val{true}, ::Val{true})::Int\n\nfirst type argument:  type of the original point (for even points -> Val{true}(), else Val{false}()) second type argument: is sum of other point indices even or odd -> (if sum is even -> Val{true}(), else Val{false}())\n\n\n\n\n\n"
},

{
    "location": "api/#SolidStateDetectors.update!-Union{Tuple{_is_weighting_potential}, Tuple{_bulk_is_ptype}, Tuple{depletion_handling_enabled}, Tuple{even_points}, Tuple{S}, Tuple{T}, Tuple{PotentialSimulationSetupRB{T,3,4,S},Int64,Val{even_points},Val{depletion_handling_enabled},Val{_bulk_is_ptype},Val{_is_weighting_potential}}} where _is_weighting_potential where _bulk_is_ptype where depletion_handling_enabled where even_points where S where T",
    "page": "API",
    "title": "SolidStateDetectors.update!",
    "category": "method",
    "text": "update!(fssrb::PotentialSimulationSetupRB{T, 3, 4, S}, RBT::DataType)::Nothing\n\nLoop over even grid points. A point is even if the sum of its cartesian indicies (of the not extended grid) is even. Even points get the red black index (rbi) = 2. ( -> rbpotential[ inds..., rbi ]).\n\n\n\n\n\n"
},

{
    "location": "api/#Documentation-1",
    "page": "API",
    "title": "Documentation",
    "category": "section",
    "text": "Modules = [SolidStateDetectors]\nOrder = [:type, :function]"
},

{
    "location": "LICENSE/#",
    "page": "LICENSE",
    "title": "LICENSE",
    "category": "page",
    "text": ""
},

{
    "location": "LICENSE/#LICENSE-1",
    "page": "LICENSE",
    "title": "LICENSE",
    "category": "section",
    "text": "using Markdown\nMarkdown.parse_file(joinpath(@__DIR__, \"..\", \"..\", \"LICENSE.md\"))"
},

]}
