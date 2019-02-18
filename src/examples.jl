
"""
    get_path_to_example_config_files()::String

Returns the path to example config files provided by the package.
"""
function get_path_to_example_config_files()::String
    return joinpath(@__DIR__, "../examples/example_config_files")
end

const SSD_examples = Dict{Symbol,String}()
SSD_examples[:Coax] = (
    joinpath(@__DIR__, "../examples/example_detector_config_files/public_Coax_config.json")
)
SSD_examples[:InvertedCoax] = (
    joinpath(@__DIR__, "../examples/example_detector_config_files/public_ivc_config.json")
)
SSD_examples[:BEGe] = (
    joinpath(@__DIR__, "../examples/example_detector_config_files/public_SegBEGe_config.json")
)
SSD_examples[:CGD] = (
    joinpath(@__DIR__, "../examples/example_detector_config_files/public_CGD_config.json")
)


