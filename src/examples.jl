
"""
    get_path_to_example_config_files()::String

Returns the path to example detector config files provided by the package.
"""
function get_path_to_example_config_files()::String
    return joinpath(@__DIR__, "../examples/example_config_files")
end

"""
    SSD_examples::Dict{Symbol,String}

containing absolute paths to example configuration files.
"""
const SSD_examples = Dict{Symbol,String}([
    :InvertedCoax   => joinpath(@__DIR__, "../examples/example_config_files/public_ivc_config.yaml"),
    :InvertedCoaxInCryostat => joinpath(@__DIR__, "../examples/example_config_files/ivc_splitted_config/public_ivc_cryostat_config.yaml"),
    :Coax           => joinpath(@__DIR__, "../examples/example_config_files/public_Coax_config.yaml"),
    :BEGe => joinpath(@__DIR__, "../examples/example_config_files/public_SegBEGe_config.yaml"),
    :CGD => joinpath(@__DIR__, "../examples/example_config_files/public_CGD_config.yaml"),
    :CGD_CylGrid => joinpath(@__DIR__, "../examples/example_config_files/public_CGD_config_cyl_grid.yaml"),
    :Spherical => joinpath(@__DIR__, "../examples/example_config_files/public_spherical_detector_config.yaml"),
    :SigGen => joinpath(@__DIR__, "../examples/example_config_files/public_ppc_config_SigGen.config"),
    :InfiniteParallelPlateCapacitor => joinpath(@__DIR__, "../examples/example_config_files/infinite_parallel_plate_capacitor.yaml"),
    :InfiniteCoaxialCapacitor => joinpath(@__DIR__, "../examples/example_config_files/infinite_coaxial_capacitor.yaml"),
    :InfiniteCoaxialCapacitorCartesianCoords => joinpath(@__DIR__, "../examples/example_config_files/infinite_coaxial_capacitor_cartesian_coords.yaml"),
    :Hexagon => joinpath(@__DIR__, "../examples/example_config_files/minimum_hexagon_config.yaml"),
    :CoaxialTorus => joinpath(@__DIR__, "../examples/example_config_files/public_coaxial_torus_config.yaml"),
])
