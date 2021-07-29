# The configuration file format was updated in v0.6.0.
# However, this configuration file still seems to be in the old format.
# 
# To update your configuration file to the new format (v0.6.0 and newer),
# open a new Julia session and load this file through
# 
#     include("<path_to_SolidStateDetectors.jl>/test/update_config_files.jl")
# 
# Afterwards, run
# 
#     update_config_file("<path_to_configuration_file>")
# 
# This method returns the file name of the updated configuration file.
# Please close the Julia session after updating the configuration files, as some
# parsing methods are overridden with old methods.

include("../src/IO/Update.jl")