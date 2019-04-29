# This file is a part of SolidStateDetectors.jl, licensed under the MIT License (MIT).

pols = Dict("neg" =>        :neg,
            "negative" =>   :neg,
            "-" =>          :neg,
            "pos" =>        :pos,
            "positive" =>   :pos,
            "+" =>          :pos)

include("Box.jl")
include("RectangularCuboid.jl")
include("Tube.jl")
include("Cone.jl")

include("plot_recipes.jl")
