struct DriftPath{T <: SSDFloat}
    e_path::Vector{<:AbstractCoordinatePoint{T}}
    h_path::Vector{<:AbstractCoordinatePoint{T}}
end

include("ChargeDriftCylindrical.jl")
include("ChargeDriftCartesian.jl")


function drift_charges( detector::SolidStateDetector{T}, grid::Grid{T, 3},
                        starting_positions::Vector{CartesianPoint{T}}, 
                        velocity_field_e::Interpolations.Extrapolation{SVector{3, Float64}, 3},
                        velocity_field_h::Interpolations.Extrapolation{SVector{3, Float64}, 3}; 
                        Δt::T = T(1f-9), n_steps::Int = 2000)::DriftPath{T} where {T <: SSDFloat}

    drift_path_e::Vector{CartesianPoint{T}} = Vector{CartesianPoint{T}}(undef, n_steps)
    drift_path_h::Vector{CartesianPoint{T}} = Vector{CartesianPoint{T}}(undef, n_steps)

    for start_pos in starting_positions
        drift_charge!(drift_path_e, detector, grid, start_pos, Δt, velocity_field_e)
        drift_charge!(drift_path_h, detector, grid, start_pos, Δt, velocity_field_h)
    end

    return DriftPath{T}( drift_path_e, drift_path_h )
end

