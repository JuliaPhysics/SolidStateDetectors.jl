# This file is a part of SolidStateDetectors.jl, licensed under the MIT License (MIT).

pols = Dict("neg" =>        :neg,
            "negative" =>   :neg,
            "-" =>          :neg,
            "pos" =>        :pos,
            "positive" =>   :pos,
            "+" =>          :pos)

include("CartesianBox3D.jl")
include("RectangularCuboid.jl")
include("Cylinder.jl")
include("Tube.jl")
include("Cone.jl")
include("SSDCone.jl")

include("plot_recipes.jl")
