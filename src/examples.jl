
"""
    get_path_to_example_config_files()::String

Returns the path to example detector config files provided by the package.
"""
function get_path_to_example_config_files()::String
    return joinpath(@__DIR__, "../examples/example_detector_config_files")
end

const SSD_examples = Dict{Symbol,String}()
SSD_examples[:Coax] = (
    joinpath(@__DIR__, "../examples/example_detector_config_files/public_Coax_config.json")
)
SSD_examples[:InvertedCoax] = (
    joinpath(@__DIR__, "../examples/example_detector_config_files/public_ivc_config.json")
)
SSD_examples[:InvertedCoaxInCryostat] = (
    joinpath(@__DIR__, "../examples/example_detector_config_files/ivc_splitted_config/public_ivc_cryostat_config.json")
)
SSD_examples[:BEGe] = (
    joinpath(@__DIR__, "../examples/example_detector_config_files/public_SegBEGe_config.json")
)
SSD_examples[:CGD] = (
    joinpath(@__DIR__, "../examples/example_detector_config_files/public_CGD_config.json")
)
SSD_examples[:Spherical] = (
    joinpath(@__DIR__, "../examples/example_detector_config_files/public_spherical_detector_config.json")
)
SSD_examples[:SigGen] = (
    joinpath(@__DIR__, "../examples/example_detector_config_files/public_ppc_config_SigGen.config")
)
SSD_examples[:InfiniteParallelPlateCapacitor] = (
    joinpath(@__DIR__, "../examples/example_detector_config_files/infinite_parallel_plate_capacitor.json")
)
SSD_examples[:InfiniteCoaxialCapacitor] = (
    joinpath(@__DIR__, "../examples/example_detector_config_files/infinite_coaxial_capacitor.json")
)
SSD_examples[:InfiniteCoaxialCapacitorCartesianCoords] = (
    joinpath(@__DIR__, "../examples/example_detector_config_files/infinite_coaxial_capacitor_cartesian_coords.json")
)
SSD_examples[:Hexagon] = (
    joinpath(@__DIR__, "../examples/example_detector_config_files/minimum_hexagon_config.json")
)
SSD_examples[:CoaxialTorus] = (
    joinpath(@__DIR__, "../examples/example_detector_config_files/public_coaxial_torus_config.json")
)
