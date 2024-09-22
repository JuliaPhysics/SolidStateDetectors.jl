
"""
    get_path_to_example_config_files()::String

Returns the path to the example detector configuration files provided by the package.

See also [`SSD_examples`](@ref).
"""
function get_path_to_example_config_files()::String
    return joinpath(dirname(@__DIR__), "examples", "example_config_files")
end

"""
    SSD_examples::Dict{Symbol,String}

Dictionary with the paths to the example detector configuration files provided by the package.

Find the possible keys of the dictionary with `keys(SSD_examples)`.

The example detector configuration files can be loaded via
```
path_to_config_file = SSD_examples[:InvertedCoax]
sim = Simulation(path_to_config_file)
```
"""
const SSD_examples = Dict{Symbol,String}([
    :TwoSpheresCapacitor                     => joinpath(get_path_to_example_config_files(), "two_spheres_capacitor.yaml"),
    :InvertedCoax                            => joinpath(get_path_to_example_config_files(), "public_ivc_config.yaml"),
    :InvertedCoaxInCryostat                  => joinpath(get_path_to_example_config_files(), "ivc_splitted_config/public_ivc_cryostat_config.yaml"),
    :InvertedCoaxTrapping                    => joinpath(get_path_to_example_config_files(), "public_ivc_trapping_config.yaml"),
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
    :IsochroneTest                           => joinpath(get_path_to_example_config_files(), "isochrone_test.yaml"),
    :ConeSym                                 => joinpath(get_path_to_example_config_files(), "Cone/cone_config_sym.yaml"),
    :Cone2D                                  => joinpath(get_path_to_example_config_files(), "Cone/cone_config_2D.yaml"),
])
