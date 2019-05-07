function clear( a::Array{T, N} )::Array{T, N} where {T, N}
    return Array{T, N}(undef, size(a) .- size(a))
end

include("RedBlack/RedBlack.jl")
include("GeometricalWeights/GeometricalWeights.jl")
include("PotentialSimulationSetups/PotentialSimulationSetups.jl")

include("ChargeDensity.jl")
include("DielectricDistribution.jl")
include("ElectricPotential.jl")
include("WeightingPotential.jl")

include("SimulationAlgorithms/Cylindrical.jl")
include("SimulationAlgorithms/Cartesian3D.jl")
include("ConvergenceAndRefinement.jl")

include("CalculateElectricPotential.jl")
include("CalculateWeightingPotential.jl")