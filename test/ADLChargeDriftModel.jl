using SolidStateDetectors: getVe, getVh, Vl, geom_sigdigits
using StaticArrays
using LinearAlgebra

config_dict = Dict(
    "High-purity germanium" => joinpath(@__DIR__,"../src/ChargeDriftModels/ADL/drift_velocity_config.json"),
    "Silicon" => joinpath(@__DIR__,"../src/ChargeDriftModels/ADL/drift_velocity_Si_300K_config.json")
)

for T in (Float32, Float64) 
    @testset "Precision type: $(T)" begin
        
        for material in keys(config_dict)
            
            @testset "$material" begin
                config = config_dict[material]
                
                @testset "x-axis aligned with <100> axis" begin

                    cdm = ADLChargeDriftModel(config, T=T, phi110 = T(π/4)); # <100> aligned with x-axis
                    E = collect(T, 10 .^(0:0.1:6));

                    r100 = @SVector T[1,0,0]
                    r111 = @SVector T[sqrt(1/3), sqrt(1/3), sqrt(1/3)]

                    @testset "Electrons <100>" begin
                        ve100_cdm = broadcast(E -> abs(getVe(E * r100, cdm) ⋅ r100), E)
                        ve100_true = map(E -> Vl(E,cdm.electrons.axis100), E)
                        @test isapprox(ve100_cdm, ve100_true, rtol = 10.0^(-geom_sigdigits(T)))
                    end

                    @testset "Electrons <111>" begin
                        ve111_cdm = broadcast(E -> abs(getVe(E * r111, cdm) ⋅ r111), E)
                        ve111_true = map(E -> Vl(E,cdm.electrons.axis111), E)
                        @test isapprox(ve111_cdm, ve111_true, rtol = 10.0^(-geom_sigdigits(T)))
                    end

                    @testset "Holes <100>" begin
                        vh100_cdm = broadcast(E -> abs(getVh(E * r100, cdm) ⋅ r100), E)
                        vh100_true = map(E -> Vl(E,cdm.holes.axis100), E)
                        @test isapprox(vh100_cdm, vh100_true, rtol = 10.0^(-geom_sigdigits(T)))
                    end

                    @testset "Holes <111>" begin
                        vh111_cdm = broadcast(E -> abs(getVh(E * r111, cdm) ⋅ r111), E)
                        vh111_true = map(E -> Vl(E,cdm.holes.axis111), E)
                        @test isapprox(vh111_cdm, vh111_true, rtol = 10.0^(-geom_sigdigits(T)))
                    end
                end
                
                @testset "x-axis and <100> axis intersect at an angle of -30°" begin

                    cdm = ADLChargeDriftModel(config, T=T, phi110 = T(π/12)); # <110> and y-axis intersect at an angle of 15°
                    E = collect(T, 10 .^(0:0.1:6));

                    r100 = @SVector T[sqrt(3/4),-1/2,0]
                    r111 = @SVector T[(1 + sqrt(1/3))/2, (1 - sqrt(1/3))/2, sqrt(1/3)]

                    @testset "Electrons <100>" begin
                        ve100_cdm = broadcast(E -> abs(getVe(E * r100, cdm) ⋅ r100), E)
                        ve100_true = map(E -> Vl(E,cdm.electrons.axis100), E)
                        @test isapprox(ve100_cdm, ve100_true, rtol = 10.0^(-geom_sigdigits(T)))
                    end

                    @testset "Electrons <111>" begin
                        ve111_cdm = broadcast(E -> abs(getVe(E * r111, cdm) ⋅ r111), E)
                        ve111_true = map(E -> Vl(E,cdm.electrons.axis111), E)
                        @test isapprox(ve111_cdm, ve111_true, rtol = 10.0^(-geom_sigdigits(T)))
                    end

                    @testset "Holes <100>" begin
                        vh100_cdm = broadcast(E -> abs(getVh(E * r100, cdm) ⋅ r100), E)
                        vh100_true = map(E -> Vl(E,cdm.holes.axis100), E)
                        @test isapprox(vh100_cdm, vh100_true, rtol = 10.0^(-geom_sigdigits(T)))
                    end

                    @testset "Holes <111>" begin
                        vh111_cdm = broadcast(E -> abs(getVh(E * r111, cdm) ⋅ r111), E)
                        vh111_true = map(E -> Vl(E,cdm.holes.axis111), E)
                        @test isapprox(vh111_cdm, vh111_true, rtol = 10.0^(-geom_sigdigits(T)))
                    end
                end
            end
        end
    
        config = joinpath(@__DIR__,"../src/ChargeDriftModels/ADL/drift_velocity_config_axes.json")

        @testset "All axes in one plane" begin

            cdm = ADLChargeDriftModel(config, T=T); # <110> and y-axis intersect at an angle of 15°
            E = collect(T, 10 .^(0:0.1:6));

            r100 = @SVector T[1,0,0]
            r111 = @SVector T[sqrt(1/3),sqrt(2/3),0]

            @testset "Electrons <100>" begin
                ve100_cdm = broadcast(E -> abs(getVe(E * r100, cdm) ⋅ r100), E)
                ve100_true = map(E -> Vl(E,cdm.electrons.axis100), E)
                @test isapprox(ve100_cdm, ve100_true, rtol = 10.0^(-geom_sigdigits(T)))
            end

            @testset "Electrons <111>" begin
                ve111_cdm = broadcast(E -> abs(getVe(E * r111, cdm) ⋅ r111), E)
                ve111_true = map(E -> Vl(E,cdm.electrons.axis111), E)
                @test isapprox(ve111_cdm, ve111_true, rtol = 10.0^(-geom_sigdigits(T)))
            end

            @testset "Holes <100>" begin
                vh100_cdm = broadcast(E -> abs(getVh(E * r100, cdm) ⋅ r100), E)
                vh100_true = map(E -> Vl(E,cdm.holes.axis100), E)
                @test isapprox(vh100_cdm, vh100_true, rtol = 10.0^(-geom_sigdigits(T)))
            end

            @testset "Holes <111>" begin
                vh111_cdm = broadcast(E -> abs(getVh(E * r111, cdm) ⋅ r111), E)
                vh111_true = map(E -> Vl(E,cdm.holes.axis111), E)
                @test isapprox(vh111_cdm, vh111_true, rtol = 10.0^(-geom_sigdigits(T)))
            end
        end
    end
end

