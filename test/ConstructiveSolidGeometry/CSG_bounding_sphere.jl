# This file is a part of SolidStateDetectors.jl, licensed under the MIT License (MIT).

using Test
using SolidStateDetectors
using LinearAlgebra

import SolidStateDetectors.ConstructiveSolidGeometry as CSG

@testset "Bounding spheres are conservative" begin
    # Point-set conservativeness: pt ∈ geometry ⇒ pt ∈ bounding sphere.
    # Surface sample points of CSG differences may lie on the (oversized)
    # subtracted solid outside the actual geometry, so filter by membership.
    T = Float32
    for (key, path) in SSD_examples
        endswith(path, ".yaml") || continue
        sim = try Simulation{T}(path) catch; continue end
        geometries = Any[sim.detector.semiconductor.geometry]
        foreach(c -> push!(geometries, c.geometry), sim.detector.contacts)
        if !ismissing(sim.detector.passives)
            foreach(p -> push!(geometries, p.geometry), sim.detector.passives)
        end
        for g in geometries
            center, radius = CSG.bounding_sphere(g)
            r2 = radius^2
            n_out = count(CSG.sample(g, T(0.003))) do p
                cp = CartesianPoint(p)
                in(cp, g) && CSG.distance_squared(cp - center) > r2
            end
            @test n_out == 0
        end
    end
end

@testset "Bounding sphere set operations" begin
    T = Float64
    box(o) = CSG.Box{T}(CSG.ClosedPrimitive; hX = T(1), hY = T(1), hZ = T(1), origin = o)
    a = box(CartesianPoint{T}(0, 0, 0))
    b = box(CartesianPoint{T}(10, 0, 0))
    ca, ra = CSG.bounding_sphere(a)
    cb, rb = CSG.bounding_sphere(b)
    cu, ru = CSG.bounding_sphere(a + b)
    # the union sphere must enclose both member spheres (tangency ⇒ ulp-scale slack)
    @test norm(cu - ca) + ra <= ru + 4 * eps(ru)
    @test norm(cu - cb) + rb <= ru + 4 * eps(ru)
    # a difference keeps the sphere of the minuend
    @test CSG.bounding_sphere(a - b) == (ca, ra)
    # an intersection may use the smaller of the two spheres
    small = CSG.Box{T}(CSG.ClosedPrimitive; hX = T(0.1), hY = T(0.1), hZ = T(0.1), origin = CartesianPoint{T}(0, 0, 0))
    ci, ri = CSG.bounding_sphere(a & small)
    @test ri <= min(ra, CSG.bounding_sphere(small)[2]) + eps(T)
end
