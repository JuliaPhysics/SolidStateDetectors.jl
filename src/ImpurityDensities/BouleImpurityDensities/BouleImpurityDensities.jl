"""
Boule impurity densities are meant to be used when the impurity density is defined in the boule coordinates, where the z-axis is aligned with the boule growth direction. 

Different models are provided. In each the field `det_z0` is the z-coordinate of the detector origin in boule coordinates. The z-direction of the detector is opposite to the z-direction of the boule coordinates.

In this matter the detector impurities are automatically determined from those of the boule, depending on `det_z0`.
"""

include("LinBouleImpurityDensity.jl")
include("LinExpBouleImpurityDensity.jl")