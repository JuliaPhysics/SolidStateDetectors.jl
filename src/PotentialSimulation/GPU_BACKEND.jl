using .CUDA
include("PotentialSimulationSetups/PotentialSimulationSetupRBGPU.jl")
include("SimulationAlgorithms/Cylindrical_GPU.jl")
include("PotentialSimulationSetups/BoundaryConditions/BoundaryConditionsCylindricalGPU.jl")