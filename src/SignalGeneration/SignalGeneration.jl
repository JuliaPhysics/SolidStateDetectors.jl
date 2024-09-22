@inline function get_interpolation(wpot::Interpolations.Extrapolation{T, 3}, pt::CartesianPoint{T}, ::Type{Cartesian})::T where {T <: SSDFloat}
    return wpot(pt.x, pt.y, pt.z)::T
end
@inline function get_interpolation(wpot::Interpolations.Extrapolation{T, 3}, pt::CartesianPoint{T}, ::Type{Cylindrical})::T where {T <: SSDFloat}
    pt_cyl::CylindricalPoint{T} = CylindricalPoint(pt)
    return wpot(pt_cyl.r, pt_cyl.φ, pt_cyl.z)::T
end
@inline function get_interpolation(wpot::Interpolations.Extrapolation{T, 3}, pt::CylindricalPoint{T}, ::Type{Cartesian})::T where {T <: SSDFloat}
    pt_car::CartesianPoint{T} = CartesianPoint(pt)
    return wpot(pt_car.x, pt_car.y, pt_car.z)::T
end
@inline function get_interpolation(wpot::Interpolations.Extrapolation{T, 3}, pt::CylindricalPoint{T}, ::Type{Cylindrical})::T where {T <: SSDFloat}
    return wpot(pt.r, pt.φ, pt.z)::T
end


function add_signal!(
        signal::AbstractVector{T}, 
        timestamps::AbstractVector{T}, 
        path::AbstractVector{CartesianPoint{T}}, 
        pathtimestamps::AbstractVector{T}, 
        charge::T,          
        wpot::Interpolations.Extrapolation{T, 3}, 
        S::CoordinateSystemType,
        ctm::AbstractChargeTrappingModel{T} = NoChargeTrappingModel{T}()
    )::Nothing where {T <: SSDFloat}

    tmp_signal::Vector{T} = _calculate_signal(ctm, path, pathtimestamps, charge, wpot, S)
    itp = interpolate!( (pathtimestamps,), tmp_signal, Gridded(Linear()))
    t_max::T = last(pathtimestamps)
    i_max::Int = searchsortedlast(timestamps, t_max)
    signal[1:i_max] += itp(timestamps[1:i_max])
    if length(signal) > i_max
        signal[i_max+1:end] .+= tmp_signal[end]
    end
    nothing
end


function add_signal!(signal::AbstractVector{T}, timestamps::AbstractVector{T}, path::EHDriftPath{T}, charge::T, 
        wpot::Interpolations.Extrapolation{T, 3}, S::CoordinateSystemType, ctm::AbstractChargeTrappingModel{T} = NoChargeTrappingModel{T}())::Nothing where {T <: SSDFloat}
    add_signal!(signal, timestamps, path.e_path, path.timestamps_e, -charge, wpot, S, ctm) # electrons induce negative charge
    add_signal!(signal, timestamps, path.h_path, path.timestamps_h,  charge, wpot, S, ctm)
    nothing
end

function add_signal!(signal::AbstractVector{T}, timestamps::AbstractVector{<:RealQuantity}, paths::Vector{<:EHDriftPath{T}}, charges::Vector{T}, 
        wpot::Interpolations.Extrapolation{T, 3}, S::CoordinateSystemType, ctm::AbstractChargeTrappingModel{T} = NoChargeTrappingModel{T}())::Nothing where {T <: SSDFloat}
    for ipath in eachindex(paths)
        add_signal!(signal, timestamps, paths[ipath], charges[ipath], wpot, S, ctm)
    end
    nothing
end
