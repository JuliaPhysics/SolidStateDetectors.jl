
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
    :InvertedCoax                            => joinpath(get_path_to_example_config_files(), "public_ivc_config.yaml"),
    :InvertedCoaxInCryostat                  => joinpath(get_path_to_example_config_files(), "ivc_splitted_config/public_ivc_cryostat_config.yaml"),
    :Coax                                    => joinpath(get_path_to_example_config_files(), "public_Coax_config.yaml"),
    :BEGe                                    => joinpath(get_path_to_example_config_files(), "public_SegBEGe_config.yaml"),
    :CGD                                     => joinpath(get_path_to_example_config_files(), "public_CGD_config.yaml"),
    :CGD_CylGrid                             => joinpath(get_path_to_example_config_files(), "public_CGD_config_cyl_grid.yaml"),
    :Spherical                               => joinpath(get_path_to_example_config_files(), "public_spherical_detector_config.yaml"),
    :SigGen                                  => joinpath(get_path_to_example_config_files(), "public_ppc_config_SigGen.config"),
    :InfiniteParallelPlateCapacitor          => joinpath(get_path_to_example_config_files(), "infinite_parallel_plate_capacitor.yaml"),
    :InfiniteCoaxialCapacitor                => joinpath(get_path_to_example_config_files(), "infinite_coaxial_capacitor.yaml"),
    :InfiniteCoaxialCapacitorCartesianCoords => joinpath(get_path_to_example_config_files(), "infinite_coaxial_capacitor_cartesian_coords.yaml"),
    :Hexagon                                 => joinpath(get_path_to_example_config_files(), "minimum_hexagon_config.yaml"),
    :CoaxialTorus                            => joinpath(get_path_to_example_config_files(), "public_coaxial_torus_config.yaml"),
])
