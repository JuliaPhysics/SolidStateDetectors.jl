struct DriftPath{T <: SSDFloat}
    e_path::Vector{<:AbstractCoordinatePoint{T}}
    h_path::Vector{<:AbstractCoordinatePoint{T}}
end

function get_velocity_vector(interpolation_field::Interpolations.Extrapolation{<:StaticVector{3}, 3}, point::CartesianPoint{T})::CartesianVector{T} where {T <: SSDFloat}  
    return CartesianVector{T}(interpolation_field(point.x, point.y, point.z))
end

@inline get_velocity_vector(interpolated_vectorfield, point::CylindricalPoint{T}) where {T <: SSDFloat} =
    CartesianVector{T}(interpolated_vectorfield(point.r, point.φ, point.z))


include("ChargeDriftCylindrical.jl")
include("ChargeDriftCartesian.jl")


function _drift_charges( detector::SolidStateDetector{T}, grid::Grid{T, 3}, point_types::PointTypes{T, 3},
                        starting_points::Vector{CartesianPoint{T}}, 
                        velocity_field_e::Interpolations.Extrapolation{SVector{3, Float64}, 3},
                        velocity_field_h::Interpolations.Extrapolation{SVector{3, Float64}, 3}; 
                        Δt::T = T(1f-9), n_steps::Int = 2000, verbose::Bool = true)::Vector{DriftPath{T}} where {T <: SSDFloat}

    drift_paths::Vector{DriftPath{T}} = Vector{DriftPath{T}}(undef, length(starting_points)) 

    for i in eachindex(starting_points)
        drift_path_e::Vector{CartesianPoint{T}} = zeros(CartesianPoint{T}, n_steps )#Vector{CartesianPoint{T}}(undef, n_steps)
        drift_path_h::Vector{CartesianPoint{T}} = zeros(CartesianPoint{T}, n_steps )#Vector{CartesianPoint{T}}(undef, n_steps)
        drift_charge!(drift_path_e, detector, point_types, grid, starting_points[i], Δt, velocity_field_e, verbose = verbose)
        drift_charge!(drift_path_h, detector, point_types, grid, starting_points[i], Δt, velocity_field_h, verbose = verbose)
        drift_paths[i] = DriftPath{T}( drift_path_e, drift_path_h )
    end

    return drift_paths 
end

function _drift_charge( detector::SolidStateDetector{T}, grid::Grid{T, 3}, point_types::PointTypes{T, 3},
                       starting_point::CartesianPoint{T}, 
                       velocity_field_e::Interpolations.Extrapolation{SVector{3, Float64}, 3},
                       velocity_field_h::Interpolations.Extrapolation{SVector{3, Float64}, 3}; 
                       Δt::T = T(1f-9), n_steps::Int = 2000, verbose::Bool = true)::Vector{DriftPath{T}} where {T <: SSDFloat}
    return _drift_charges(detector, grid, point_types, [starting_point], velocity_field_e, velocity_field_h, Δt = Δt, n_steps = n_steps, verbose = verbose)
end

