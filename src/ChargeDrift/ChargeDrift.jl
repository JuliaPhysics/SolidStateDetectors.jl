struct DriftPath{T <: SSDFloat}
    e_path::Vector{<:AbstractCoordinatePoint{T}}
    h_path::Vector{<:AbstractCoordinatePoint{T}}
end

include("ChargeDriftCylindrical.jl")
include("ChargeDriftCartesian.jl")


function drift_charges( detector::SolidStateDetector{T}, grid::Grid{T, 3},
                        starting_points::Vector{CartesianPoint{T}}, 
                        velocity_field_e::Interpolations.Extrapolation{SVector{3, Float64}, 3},
                        velocity_field_h::Interpolations.Extrapolation{SVector{3, Float64}, 3}; 
                        Δt::T = T(1f-9), n_steps::Int = 2000)::Vector{DriftPath{T}} where {T <: SSDFloat}

    drift_paths::Vector{DriftPath{T}} = Vector{DriftPath{T}}(undef, length(starting_points))
    drift_path_e::Vector{CartesianPoint{T}} = Vector{CartesianPoint{T}}(undef, n_steps)
    drift_path_h::Vector{CartesianPoint{T}} = Vector{CartesianPoint{T}}(undef, n_steps)

    for i in eachindex(starting_points)
        drift_charge!(drift_path_e, detector, grid, starting_points[i], Δt, velocity_field_e)
        drift_charge!(drift_path_h, detector, grid, starting_points[i], Δt, velocity_field_h)
        drift_paths[i] = DriftPath{T}( drift_path_e, drift_path_h )
    end

    return drift_paths 
end

function drift_charge( detector::SolidStateDetector{T}, grid::Grid{T, 3},
                       starting_point::CartesianPoint{T}, 
                       velocity_field_e::Interpolations.Extrapolation{SVector{3, Float64}, 3},
                       velocity_field_h::Interpolations.Extrapolation{SVector{3, Float64}, 3}; 
                       Δt::T = T(1f-9), n_steps::Int = 2000)::Vector{DriftPath{T}} where {T <: SSDFloat}

    drift_path_e::Vector{CartesianPoint{T}} = Vector{CartesianPoint{T}}(undef, n_steps)
    drift_path_h::Vector{CartesianPoint{T}} = Vector{CartesianPoint{T}}(undef, n_steps)
    drift_charge!(drift_path_e, detector, grid, starting_point, Δt, velocity_field_e)
    drift_charge!(drift_path_h, detector, grid, starting_point, Δt, velocity_field_h)
    return DriftPath{T}( drift_path_e, drift_path_h )
end

