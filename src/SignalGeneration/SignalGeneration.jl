@inline function get_interpolation(wp::Interpolations.Extrapolation{T, 3}, pt::CartesianPoint{T}, S::Val{:cartesian})::T where {T <: SSDFloat}
    return wp(pt.x, pt.y, pt.z)::T
end
@inline function get_interpolation(wp::Interpolations.Extrapolation{T, 3}, pt::CartesianPoint{T}, S::Val{:cylindrical})::T where {T <: SSDFloat}
    p::CylindricalPoint{T} = CylindricalPoint(pt)
    return wp(p.r, p.φ, p.z)::T
end


function add_signal!(signal::Vector{T}, path::Vector{CartesianPoint{T}}, charge::T, wp::Interpolations.Extrapolation{T, 3}, S::Union{Val{:cylindrical}, Val{:cartesian}})::Nothing where {T <: SSDFloat}
    @inbounds for i in eachindex(signal)
        signal[i] += get_interpolation(wp, path[i], S)::T * charge
    end
    nothing
end
function add_signal!(signal::Vector{T}, path::DriftPath{T}, charge::T, wp::Interpolations.Extrapolation{T, 3}, S::Union{Val{:cylindrical}, Val{:cartesian}})::Nothing where {T <: SSDFloat}
    add_signal!(signal, path.e_path, -charge, wp, S) # electrons induce negative charge
    add_signal!(signal, path.h_path,  charge, wp, S)
    nothing
end
function add_signal!(signal::Vector{T}, paths::Vector{DriftPath{T}}, charges::Vector{T}, wp::Interpolations.Extrapolation{T, 3}, S::Union{Val{:cylindrical}, Val{:cartesian}})::Nothing where {T <: SSDFloat}
    for ipath in eachindex(paths)
        add_signal!(signal, paths[ipath], charges[ipath], wp, S)
    end
    nothing
end

function signal_contributions_from_drift_paths(drift_paths, energy_depositions::AbstractVector{T}, Wpot_interp::Interpolations.Extrapolation{T,3}) where T<:Real
    charge_signal_e = zeros(T,size(drift_paths[1].e_path,1))
    charge_signal_h = zeros(T,size(drift_paths[1].e_path,1))
    @inbounds for i in eachindex(drift_paths[1].e_path)
        for (idp, drift_path) in enumerate(drift_paths)
            # charge_signal_e[i] = muladd(energy_depositions[idp] ,charge_signal_e[i], _get_wpot_at(Wpot_interp, drift_path.e_path[i]))
            # charge_signal_h[i] = muladd(energy_depositions[idp] ,charge_signal_h[i], _get_wpot_at(Wpot_interp, drift_path.h_path[i]))
            charge_signal_e[i] +=  _get_wpot_at(Wpot_interp, drift_path.e_path[i]) * energy_depositions[idp]
            charge_signal_h[i] +=  _get_wpot_at(Wpot_interp, drift_path.h_path[i]) * energy_depositions[idp]
        end
    end
    return charge_signal_e .* -1, charge_signal_h
end
@deprecate signal_contributions_from_drift_paths(drift_paths, energy_depositions::AbstractVector{T}, Wpot_interp::Interpolations.Extrapolation{T,3}) where {T <: SSDFloat} get_signal(dp::DriftPath{T}, energy::T, wp::Interpolations.Extrapolation{T,3,ITPT,IT,ET} where ET where IT where ITPT) where {T <: SSDFloat}