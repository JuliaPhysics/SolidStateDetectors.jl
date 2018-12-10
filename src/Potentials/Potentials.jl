const electron_charge = Float64(1.602176487e-19)
const epsilon_0 = Float64(8.8541878176e-12)
const epsilon_Ge = Float64(16)

const effective_electron_charge_germanium = Float64(electron_charge / (epsilon_0 * epsilon_Ge))
const effective_electron_charge_vacuum = Float64(electron_charge / epsilon_0 )

include("PrecalculatedWeights.jl")

include("RedBlackCylindricalAlgorithm.jl")

include("ConvergenceAndRefinement.jl")

include("CalculateElectricPotential.jl")
include("CalculateWeightingPotential.jl")