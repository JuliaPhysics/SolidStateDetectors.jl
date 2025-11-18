include("PowerLawTemperatureModel.jl")
include("VacuumTemperatureModel.jl")

scale_to_temperature(cdm::AbstractChargeDriftModel, ::Missing) = cdm