include("LinearModel.jl")
include("BoltzmannModel.jl")
include("PowerLawModel.jl")
include("VacuumModel.jl")
include("SquareRootModel.jl")

show(io::IO, tm::AbstractTemperatureModel) = print(io, tm) 
show(io::IO,::MIME"text/plain", tm::AbstractTemperatureModel) = show(io, tm)

scale_to_given_temperature(Emag::T, m::AbstractTemperatureModel{T}) where T <: SSDFloat = scale_to_given_temperature(m)