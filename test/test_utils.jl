# This file is a part of SolidStateDetectors.jl, licensed under the MIT License (MIT).

using TimerOutputs
# using CUDAKernels # Uncomment this line in order to run all tests on (CUDA) GPU

device_array_type = (@isdefined CUDAKernels) ? CUDAKernels.CUDA.CuArray : Array

testtimer() = get_timer("_default_testtimer_")

macro timed_testset(title, body)
    quote
        tmr = testtimer()
        _title = $(esc(title))
        @timeit tmr "$_title" begin
            @testset "$_title" begin
                $(esc(body))
            end
        end
    end
end


function timed_calculate_electric_potential!(args...; kwargs...)
    @timeit testtimer() "calculate_electric_potential!" begin
        calculate_electric_potential!(args...; kwargs...)
    end
end

function timed_calculate_electric_field!(args...; kwargs...)
    @timeit testtimer() "calculate_electric_field!" begin
        calculate_electric_field!(args...; kwargs...)
    end
end

function timed_calculate_weighting_potential!(args...; kwargs...)
    @timeit testtimer() "calculate_weighting_potential!" begin
        calculate_weighting_potential!(args...; kwargs...)
    end
end

function timed_simulate_waveforms(args...; kwargs...)
    @timeit testtimer() "simulate_waveforms" begin
        simulate_waveforms(args...; kwargs...)
    end
end

function timed_simulate!(args...; kwargs...)
    @timeit testtimer() "simulate!" begin
        simulate!(args...; kwargs...)
    end
end

function timed_calculate_capacitance_matrix(args...; kwargs...)
    @timeit testtimer() "calculate_capacitance_matrix" begin
        calculate_capacitance_matrix(args...; kwargs...)
    end
end

function timed_estimate_depletion_voltage(args...; kwargs...)
    @timeit testtimer() "estimate_depletion_voltage" begin
        estimate_depletion_voltage(args...; kwargs...)
    end
end
