# This file is a part of SolidStateDetectors.jl, licensed under the MIT License (MIT).

using Test
using SolidStateDetectors

@info "Running tests with $(Base.Threads.nthreads()) Julia threads active."

include("test_utils.jl")

@timed_testset "Comparison to analytic solutions" begin
    include("test_analytic_solutions.jl")
end

@timed_testset "SOR GPU Backend" begin
    include("test_gpu_backend.jl")
end

@timed_testset "ConstructiveSolidGeometry" begin
    include("ConstructiveSolidGeometry/CSG_IO.jl")
    include("ConstructiveSolidGeometry/CSG_primitives.jl")
end

@timed_testset "Test real detectors" begin
    include("test_real_detectors.jl")
end

@timed_testset "Test charge drift" begin
    include("test_charge_drift_models.jl")
end

@timed_testset "Fano factor" begin
    

end

@timed_testset "Isochrone" begin
    include("test_isochrone.jl")
end

@timed_testset "IO" begin
    include("test_io.jl")
end 

@timed_testset "Geant4 extension" begin
    if Sys.WORD_SIZE == 64 include("test_geant4.jl") end
end

@timed_testset "Depletion estimation" begin
    include("test_depletion.jl")
end

@timed_testset "Two step potential refinement" begin
    include("test_refine_existing_potential.jl")
end

display(testtimer())