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
    "text": "Example minimum config file for an Inverted Coax detector (IVC) plus explanations.  Remember, comments are not allowed in JSON files and have to be deleted if you want to use it.{\n    \"name\":\"ExampleInvertedCoax\", // Arbitrary name of the detector\n    \"class\":\"InvertedCoax\",  // either \"Coax\", \"BEGe\", \"InvertedCoax\"\n    \"type\":\"p\", // either \"p\", \"ptype\", \"p-type\", \"n\", \"ntype\" or \"n-type\"\n    \"cyclic\":0, // The periodicity of the detector in degree. \n                // `0` means complete symmetric in ϕ -> 2D simulation. \n                // The condition \"360 / `cyclic` = n; n = 1, 2, 3, ...\" must be fulfilled.\n    \"materials\":{\n        \"detector\":\"HPGe\",  // Material of the detector \n        \"environment\":\"Vacuum\" // Material of the environment\n    },\n    \"geometry\":\n    {\n        \"unit\": \"mm\", // unit of the values which the user specifies in this config file\n        \"crystal\":\n        {\n            \"length\": 80.0, // crystal goes from z=0mm to z=80mm\n            \"radius\": 35.0  // crystal goes from r=0mm to r=35mm\n        },\n        \"point_contact\":\n        {\n            \"endplate\":\"bot\",\n            \"depth\":3e-4,\n            \"radius\":9.5\n        },\n        \"groove\":\n        {\n            \"endplate\":\"bot\",\n            \"rInner\":13.5,\n            \"width\":0.0,\n            \"depth\":0.0\n        },\n        \"borehole\": // IVC: Borehole starts at the top of the crystal \n        {\n          \"length\":55.0,\n          \"radius\":5.0\n        },\n        \"taper\": \n        {\n            \"outer\":\n            {\n                \"rInner\":24.42,\n                \"length\":60\n            },\n            \"inner\":\n            {\n                \"rOuter\":0,\n                \"length\":0\n            }\n        }\n    },\n    \"charge_carrier_density\":\n    {\n        \"unit\": \"cm-3\",\n        \"top\": 0.94, // times 10^10\n        \"bot\": 1.10  // just the absolute value. The sign is determined through the bulk material\n    },\n    \"segmentation\":\n    {\n        \"n_contacts_total\": 3, // n-contacts + p-contacts\n        \"core\": // core = n-contact of the detector\n        {\n            \"type\": \"Tubs\", // \"Box\" volume in cylindrical coordinates\n            \"rStart\": 0.0,\n            \"rStop\": 3.0,\n            \"phiStart\": 0,\n            \"phiStop\": 360,\n            \"zStart\": 0.0,\n            \"zStop\": 0.0,\n            \"potential\": 0.0 // Bias voltage of the core\n        },        \n        \"n_individual_segments\": 2, // now the p-contacts\n        \"S1\":\n        {\n            \"type\": \"Tubs\",\n            \"rStart\": 15.0,\n            \"rStop\": 35.0,\n            \"phiStart\": 0,\n            \"phiStop\": 360,\n            \"zStart\": 0,\n            \"zStop\": 0,\n            \"potential\": 3500.0, // Bias voltage of this segment\n            \"boundaryWidth\":\n            {\n                \"radial\":0.0,\n                \"vertical\": 0.0,\n                \"horizontal\": 0.0\n            },\n            \"repetitive\":false,\n            \"repetitions\":\n            {\n                \"radial\":0,\n                \"vertical\": 0,\n                \"horizontal\": 5\n            }\n        },\n        \"S2\":\n        {\n            \"type\": \"Tubs\",\n            \"rStart\": 35.0,\n            \"rStop\": 35.0,\n            \"phiStart\": 0,\n            \"phiStop\": 360,\n            \"zStart\": 0.0,\n            \"zStop\": 20.0,\n            \"potential\": 3500.0, // Bias voltage of this segment\n            \"boundaryWidth\":\n            {\n                \"radial\":0.0,\n                \"vertical\": 0.0,\n                \"horizontal\": 0.0\n            },\n            \"repetitive\":false,\n            \"repetitions\":\n            {\n                \"radial\": 0,\n                \"vertical\": 0,\n                \"horizontal\": 5\n            }\n        },\n        \"custom_grouping\":false, // Multiple geometrical segments can build up one contact:\n                                 // Specify which geometrical segments belong to which contact\n                                 // To specify geometrical segments `3`, `4`, `5`, `6` use \"3-4-5-6\" (not only \"3-6\")\n        \"Chn1\":\"1\",    // contact \"Chn1\" is build through the geometrical segment `1`. (This is the core in this case)\n        \"Chn2\":\"2-3\"   // contact \"Chn2\" is build through the geometrical segments `2` & `3`. \n    }\n}"
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
    "text": "To simulate a detector the detector must first be defined (Detectors). using Plots; pyplot()\nusing SolidStateDetectors\n\ndetector = SolidStateDetector(SSD_examples[:InvertedCoax]);Then, the electric potential can be calculated through the function calculate_electric_potential. It returns the electric potential E_pot in form of the the struct SolidStateDetectors.CylindricalGrid and an 3-dimensional array of point types point_types.E_pot, point_types = calculate_electric_potential(detector);A plot recipe for the struct SolidStateDetectors.CylindricalGrid exists so the result can be visualized throughplot(E_pot, φ=0)Since this is a fully φ-symmetric detector, the 3-dimensional grid has only 1 tick in the polar coordinate φ. But due to effects of the crystal axes the grid has to be extended in φ in order to deterime the electric  and drift fields. This can be done through SolidStateDetectors.extent_2D_grid_to_3D!:SolidStateDetectors.extent_2D_grid_to_3D!(E_pot, 36);Now it can also be plotted in the r-φ plane: plot(E_pot, z=0.04)"
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
    "text": "Right now (electric & weighting) potentials and point types can easily be saved and loaded via HDF5.jl.using SolidStateDetectors, HDF5\ndetector = SolidStateDetector(SSD_examples[:InvertedCoax])\nE_pot, point_types = calculate_electric_potential(detector)\n\n# write to HDF5\nh5f = h5open(\"InvertedCoaxSimulation.hdf5\", \"w\") \ng_E_pot = g_create(h5f, \"Electric Potential\")\ng_point_types = g_create(h5f, \"Point Types\")\nSolidStateDetectors.write_to_hdf5(g_E_pot, E_pot) \nSolidStateDetectors.write_to_hdf5(g_point_types, point_types) \nclose(h5f)using SolidStateDetectors, HDF5\n\n# read from HDF5\nh5f = h5open(\"InvertedCoaxSimulation.hdf5\", \"r\") \ng_E_pot = g_open(h5f, \"Electric Potential\")\ng_point_types = g_open(h5f, \"Point Types\")\nE_pot = SolidStateDetectors.read_from_hdf5(g_E_pot, ElectricPotential)\npoint_types = SolidStateDetectors.read_from_hdf5(g_point_types, PointTypes)\nclose(h5f)"
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
    "location": "api/#SolidStateDetectors.CylindricalGrid",
    "page": "API",
    "title": "SolidStateDetectors.CylindricalGrid",
    "category": "type",
    "text": "CylindricalGrid{T<:AbstractFloat} <: Grid\n\nStores the electric potential on a three dimensional cylindrical grid and also the position of the grid points.\n\nFields\n\ncyclic::T: Stores the periodicity of the system. E.g. 2π.\nr::T: The r-axis of the grid in m.\nφ::T: The φ-axis of the grid in rad.\nz::T: The r-axis of the grid in m.\npotential::Array{T, 3}: The potential values in V of each grid point.\n\n\n\n\n\n"
},

{
    "location": "api/#SolidStateDetectors.PointTypes",
    "page": "API",
    "title": "SolidStateDetectors.PointTypes",
    "category": "type",
    "text": "PointTypes{T <: AbstractFloat} <: AbstractPointTypes\n\nStores the point type information in a three dimensional cylindrical grid and also the position of the grid points.\n\nFields\n\nr::T: The r-axis of the grid in m.\nφ::T: The φ-axis of the grid in rad.\nz::T: The r-axis of the grid in m.\npointtypes::Array{PointType, 3}: The point type of each grid point.\n\n\n\n\n\n"
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
    "location": "api/#SolidStateDetectors.add_fano_noise",
    "page": "API",
    "title": "SolidStateDetectors.add_fano_noise",
    "category": "function",
    "text": "add_fano_noise(E_dep::RealQuantity, E_ionisation::RealQuantity, f_fano::Real)::RealQuantity\n\nAdd Fano noise to an energy deposition E_dep, assuming a detector material ionisation energy E_ionisation and a Fano factor f_fano.\n\n\n\n\n\n"
},

{
    "location": "api/#SolidStateDetectors.calculate_electric_potential-Tuple{SolidStateDetector}",
    "page": "API",
    "title": "SolidStateDetectors.calculate_electric_potential",
    "category": "method",
    "text": "calculate_electric_potential(det::SolidStateDetector; <keyword arguments>)::Tuple{<:Grid, PointTypes}\n\nCompute the electric potential for the given Detector det on an adaptive grid through successive over relaxation.\n\nThere are serveral <keyword arguments> which can be used to tune the computation:\n\nKeywords\n\ncoordinates::Symbol: the kind of the coordinate system of the grid. Right now only :cylindrical is possible.\nconvergence_limit::Real: convergence_limit times the bias voltage sets the convergence limit of the relaxation. The convergence value is the absolute maximum difference of the potential between two iterations of all grid points. Default of convergence_limit is 5e-6 (times bias voltage).\nmax_refinements::Int: number of maximum refinements. Default is 2. Set it to 0 to switch off refinement.\nrefinement_limits::Vector{Real}: vector of refinement limits for each dimension (in case of cylindrical coordinates the order is r, φ, z). A refinement limit (e.g. refinement_limits[1]) times the bias voltage of the detector det is the maximum allowed voltage difference between two neighbouring grid points in the respective dimension. When the difference is larger, new points are created inbetween. Default is [1e-4, 1e-4, 1e-4].\nmin_grid_spacing::Vector{Real}: vector of the mimimum allowed distance between two grid points for each dimension. For normal coordinates the unit is meter. For angular coordinates, the unit is radiance. It prevents the refinement to make the grid to fine. Default is [1e-4, 1e-2, 1e-4].\ndepletion_handling::Bool: enables the handling of undepleted regions. Default is false.\nnthreads::Int: Number of threads to use in the computation. Default is Base.Threads.nthreads(). The environment variable JULIA_NUM_THREADS must be set appropriately before the Julia session was started (e.g. export JULIA_NUM_THREADS=8 in case of bash).\nsor_consts::Vector{<:Real}: Two element array. First element contains the SOR constant for r = 0. Second contains the constant at the outer most grid point in r. A linear scaling is applied in between. First element should be smaller than the second one and both should be ∈ [1.0, 2.0]. Default is [1.4, 1.85].\nreturn_2π_grid::Bool: If set to true, the grid is extended to 2π in φ after it is computated. Default is true. This keyword is ignored if the simulation is in 2D. Use extent_2D_grid_to_3D() function to extend a 2D grid into 3D.\nmax_n_iterations::Int: Set the maximum number of iterations which are performed after each grid refinement. Default is 10000. If set to -1 there will be no limit.\nverbose::Bool=true: Boolean whether info output is produced or not.\ninit_grid_spacing::Vector{<:Real}: Initial spacing of the grid. Default is [2e-3, 2π / 72, 2e-3] <=> [2mm, 5 degree, 2mm ]\n\nAdditional Information\n\nThe function always turns the detector into a p-type detector to compute the potential. At the end it turns the signs to make it n-type again if it was n-type.\n\n\n\n\n\n"
},

{
    "location": "api/#SolidStateDetectors.calculate_weighting_potential-Tuple{SolidStateDetector,Int64}",
    "page": "API",
    "title": "SolidStateDetectors.calculate_weighting_potential",
    "category": "method",
    "text": "calculate_weighting_potential(det::SolidStateDetector, channels::Array{Int, 1}; <keyword arguments>)::Grid\n\nCompute the weighting potential of the channels of the given Detector det on an adaptive grid through successive over relaxation. channels is a list of the channels which are fixed to 1. All other channels are fixed to 0.\n\nThere are serveral <keyword arguments> which can be used to tune the computation:\n\nKeywords\n\ncoordinates::Symbol: the kind of the coordinate system of the grid. Right now only :cylindrical is possible.\nconvergence_limit::Real: convergence_limit times the bias voltage sets the convergence limit of the relaxation. The convergence value is the absolute maximum difference of the potential between two iterations of all grid points. Default of convergence_limit is 5e-6 (times bias voltage).\nmax_refinements::Int: number of maximum refinements. Default is 2. Set it to 0 to switch off refinement.\nrefinement_limits::Vector{Real}: vector of refinement limits for each dimension (in case of cylindrical coordinates the order is r, φ, z). A refinement limit (e.g. refinement_limits[1]) times the bias voltage of the detector det is the maximum allowed voltage difference between two neighbouring grid points in the respective dimension. When the difference is larger, new points are created inbetween. Default is [1e-4, 1e-4, 1e-4].\nmin_grid_spacing::Vector{Real}: vector of the mimimum allowed distance between two grid points for each dimension. For normal coordinates the unit is meter. For angular coordinates the unit is radiance. It prevents the refinement to make the grid to fine. Default is [1e-4, 1e-2, 1e-4].\ndepletion_handling::Bool: NOT IMPLEMENTET YET. Enables the handling of undepleted regions. Default is false.\nnthreads::Int: Number of threads to use in the computation. Default is Base.Threads.nthreads(). The environment variable JULIA_NUM_THREADS must be set appropriately before the Julia session was started (e.g. export JULIA_NUM_THREADS=8 in case of bash).\nsor_consts::Vector{<:Real}: Two element array. First element contains the SOR constant for r = 0. Second contains the constant at the outer most grid point in r. A linear scaling is applied in between. First element should be smaller than the second one and both should be ∈ [1.0, 2.0].\nmax_n_iterations::Int: Set the maximum number of iterations which are performed after each grid refinement. Default is 10000. If set to -1 there will be no limit.\nverbose::Bool=true: Boolean whether info output is produced or not.\ninit_grid_spacing::Vector{<:Real}: Initial spacing of the grid. Default is [2e-3, 2π / 72, 2e-3] <=> [2mm, 5 degree, 2mm ]\n\nAdditional Information\n\n\n\n\n\n"
},

{
    "location": "api/#SolidStateDetectors.get_active_volume-Union{Tuple{T}, Tuple{CylindricalGrid,PointTypes{T}}} where T",
    "page": "API",
    "title": "SolidStateDetectors.get_active_volume",
    "category": "method",
    "text": "get_active_volume(grid::CylindricalGrid, pts::PointTypes{T}) where {T}\n\nReturns an approximation of the active volume of the detector by summing up the cell volumes of all depleted cells.\n\n\n\n\n\n"
},

{
    "location": "api/#SolidStateDetectors.PointType",
    "page": "API",
    "title": "SolidStateDetectors.PointType",
    "category": "type",
    "text": "const PointType = UInt8\n\nStores certain information about a grid point via bit-flags. \n\nRight now there are:\n\n`const update_bit      = 0x01`\n`const undepleted_bit  = 0x02`\n`const bubble_bit      = 0x04`\n`const pn_junction_bit = 0x08`\n\nHow to get information out of a PointType variable pt:\n\npt & update_bit == 0 -> do not update this point (for fixed points)     \npt & update_bit >  0 -> do update this point    \npt & undepleted_bit > 0 -> this point is undepleted\npt & bubble_bit > 0     -> this point is in an \"bubble\"-area/volume. So an undepleted region which is not in contact with an electrode (fixed points).\npt & pn_junction_bit > 0  -> this point to the solid state detector. So it is in the volume of the n-type or p-type material.\n\n\n\n\n\n"
},

{
    "location": "api/#SolidStateDetectors.extent_2D_grid_to_3D!-Tuple{CylindricalGrid,Int64}",
    "page": "API",
    "title": "SolidStateDetectors.extent_2D_grid_to_3D!",
    "category": "method",
    "text": "extent_2D_grid_to_3D!(grid::CylindricalGrid, n::Int)::Nothing\n\nThis function extend a 2-dimensional grid (only one tick in φ) to an 3-dimensional grid with n ticks in φ (symmetrically distributed up to 2π).\n\n\n\n\n\n"
},

{
    "location": "api/#SolidStateDetectors.extent_2D_grid_to_3D-Tuple{CylindricalGrid,Int64}",
    "page": "API",
    "title": "SolidStateDetectors.extent_2D_grid_to_3D",
    "category": "method",
    "text": "extent_2D_grid_to_3D(grid::CylindricalGrid, n::Int)::CylindricalGrid\n\nThis function returns a extended grid of a 2-dimensional grid grid. The extended grid has n ticks in φ (symmetrically distributed up to 2π).\n\n\n\n\n\n"
},

{
    "location": "api/#SolidStateDetectors.get_path_to_example_config_files-Tuple{}",
    "page": "API",
    "title": "SolidStateDetectors.get_path_to_example_config_files",
    "category": "method",
    "text": "get_path_to_example_config_files()::String\n\nReturns the path to example config files provided by the package.\n\n\n\n\n\n"
},

{
    "location": "api/#SolidStateDetectors.read_from_hdf5-Tuple{Any,Type{CylindricalGrid}}",
    "page": "API",
    "title": "SolidStateDetectors.read_from_hdf5",
    "category": "method",
    "text": "read_from_hdf5(input, ::Type{ElectricPotential})::ElectricPotential\n\ninput should be and HDF5.HDF5Group.\n\n\n\n\n\n"
},

{
    "location": "api/#SolidStateDetectors.read_from_hdf5-Tuple{Any,Type{PointTypes}}",
    "page": "API",
    "title": "SolidStateDetectors.read_from_hdf5",
    "category": "method",
    "text": "read_from_hdf5(input, ::Type{PointTypes})::PointTypes\n\ninput should be and HDF5.HDF5Group.\n\n\n\n\n\n"
},

{
    "location": "api/#SolidStateDetectors.write_to_hdf5-Tuple{Any,CylindricalGrid}",
    "page": "API",
    "title": "SolidStateDetectors.write_to_hdf5",
    "category": "method",
    "text": "write_to_hdf5(output, ep::ElectricPotential)\n\noutput should be and HDF5.HDF5Group.\n\n\n\n\n\n"
},

{
    "location": "api/#SolidStateDetectors.write_to_hdf5-Tuple{Any,PointTypes}",
    "page": "API",
    "title": "SolidStateDetectors.write_to_hdf5",
    "category": "method",
    "text": "write_to_hdf5(output, ep::PointTypes)\n\noutput should be and HDF5.HDF5Group.\n\n\n\n\n\n"
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
