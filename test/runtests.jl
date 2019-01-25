# This file is a part of SolidStateDetectors.jl, licensed under the MIT License (MIT).
using SolidStateDetectors
using Test

@testset "Package SolidStateDetectors" begin

    @testset "Load Config Files" begin
        @test begin
            global det_coax
            det_coax = SolidStateDetector(SSD_examples[:Coax])
            true
        end
        @test begin
            global det_bege
            det_bege = SolidStateDetector(SSD_examples[:BEGe])
            true
        end
        @test begin 
            global det_ivc
            det_ivc = SolidStateDetector(SSD_examples[:InvertedCoax])
            true
        end
    end 

    @testset "Electric Potential Simulation" begin
        @test begin
            global eps_coax
            eps_coax = calculate_electric_potential(det_coax, max_refinements = 1, convergence_limit=5e-5);
            !(any(isnan, eps_coax.potential) || any(isinf, eps_coax.potential)) 
        end
        @test begin
            global eps_bege
            eps_bege = calculate_electric_potential(det_bege, max_refinements = 1, convergence_limit=5e-5);
            !(any(isnan, eps_bege.potential) || any(isinf, eps_bege.potential)) 
        end
        @test begin
            global eps_ivc
            eps_ivc = calculate_electric_potential(det_ivc, max_refinements = 1, convergence_limit=5e-5);
            !(any(isnan, eps_ivc.potential) || any(isinf, eps_ivc.potential)) 
        end
    end
    @testset "Weighting Potential Simulation" begin
        @test begin
            global wps_coax
            wps_coax = calculate_weighting_potential(det_coax, 1, max_refinements = 1, convergence_limit=5e-5);
            !(any(isnan, wps_coax.potential) || any(isinf, wps_coax.potential))
        end
        @test begin
            global wps_bege
            wps_bege = calculate_weighting_potential(det_bege, 1, max_refinements = 1, convergence_limit=5e-5);
            !(any(isnan, wps_bege.potential) || any(isinf, wps_bege.potential)) 
        end
        @test begin
            global wps_ivc
            wps_ivc = calculate_weighting_potential(det_ivc, 1, max_refinements = 1, convergence_limit=5e-5);
            !(any(isnan, wps_ivc.potential) || any(isinf, wps_ivc.potential)) 
        end
    end

end