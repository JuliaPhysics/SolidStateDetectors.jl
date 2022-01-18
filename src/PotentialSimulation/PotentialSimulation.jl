function clear( a::Array{T, N} )::Array{T, N} where {T, N}
    return Array{T, N}(undef, size(a) .- size(a))
end

include("RedBlack/RedBlack.jl")
include("GeometricalWeights/GeometricalWeights.jl")
include("Painting/Painting.jl")
include("PotentialSimulationSetups/PotentialSimulationSetups.jl")

include("EffectiveChargeDensity.jl")
include("DielectricDistribution.jl")
include("ElectricPotential.jl")
include("WeightingPotential.jl")
include("ScalarPotential.jl")

include("SimulationAlgorithms/SOR.jl")
include("SimulationAlgorithms/CPU_Cylindrical.jl")
include("SimulationAlgorithms/CPU_Cartesian3D.jl")
include("SimulationAlgorithms/GPU_Cylindrical.jl")
include("SimulationAlgorithms/GPU_Cartesian3D.jl")
include("ConvergenceAndRefinement.jl")