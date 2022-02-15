function clear( a::Array{T, N} )::Array{T, N} where {T, N}
    return Array{T, N}(undef, size(a) .- size(a))
end

include("RedBlack.jl")
include("GeometricalWeights.jl")
include("Painting.jl")
include("PotentialSimulationSetups/PotentialSimulationSetups.jl")

include("EffectiveChargeDensity.jl")
include("DielectricDistribution.jl")
include("ElectricPotential.jl")
include("WeightingPotential.jl")
include("ScalarPotential.jl")

include("SuccessiveOverRelaxation/SOR.jl")
include("Refinement.jl")