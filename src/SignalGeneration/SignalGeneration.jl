@inline function get_interpolation(wp::Interpolations.Extrapolation{T, 3}, pt::CartesianPoint{T}, ::Type{Cartesian})::T where {T <: SSDFloat}
    return wp(pt.x, pt.y, pt.z)::T
end
@inline function get_interpolation(wp::Interpolations.Extrapolation{T, 3}, pt::CartesianPoint{T}, ::Type{Cylindrical})::T where {T <: SSDFloat}
    p::CylindricalPoint{T} = CylindricalPoint(pt)
    return wp(p.r, p.φ, p.z)::T
end
@inline function get_interpolation(wp::Interpolations.Extrapolation{T, 3}, pt::CylindricalPoint{T}, ::Type{Cartesian})::T where {T <: SSDFloat}
    p::CartesianPoint{T} = CartesianPoint(pt)
    return wp(p.x, p.y, p.z)::T
end
@inline function get_interpolation(wp::Interpolations.Extrapolation{T, 3}, pt::CylindricalPoint{T}, ::Type{Cylindrical})::T where {T <: SSDFloat}
    return wp(pt.r, pt.φ, pt.z)::T
end


function add_signal!(signal::AbstractVector{T}, timestamps::AbstractVector{T}, path::Vector{CartesianPoint{T}}, pathtimestamps::AbstractVector{T}, charge::T, 
                        wp::Interpolations.Extrapolation{T, 3}, S::CoordinateSystemType)::Nothing where {T <: SSDFloat}
    tmp_signal = Vector{T}(undef, length(pathtimestamps))
    @inbounds for i in eachindex(tmp_signal)
        tmp_signal[i] = get_interpolation(wp, path[i], S)::T * charge
    end
    itp = interpolate( (pathtimestamps,), tmp_signal, Gridded(Linear()))
    t_max::T = last(pathtimestamps)
    i_max::Int = searchsortedlast(timestamps, t_max)
    signal[1:i_max] += itp(timestamps[1:i_max])
    if length(signal) > i_max
        signal[i_max+1:end] .+= tmp_signal[end]
    end
    nothing
end


function add_signal!(signal::AbstractVector{T}, timestamps::AbstractVector{TT}, path::EHDriftPath{T, TT}, charge::T, wp::Interpolations.Extrapolation{T, 3}, S::CoordinateSystemType)::Nothing where {T <: SSDFloat, TT}
    add_signal!(signal, timestamps, path.e_path, path.timestamps_e, -charge, wp, S) # electrons induce negative charge
    add_signal!(signal, timestamps, path.h_path, path.timestamps_h,  charge, wp, S)
    nothing
end

function add_signal!(signal::AbstractVector{T}, timestamps::AbstractVector{<:RealQuantity}, paths::Vector{<:EHDriftPath{T}}, charges::Vector{T}, wp::Interpolations.Extrapolation{T, 3}, S::CoordinateSystemType)::Nothing where {T <: SSDFloat}
    for ipath in eachindex(paths)
        add_signal!(signal, timestamps, paths[ipath], charges[ipath], wp, S)
    end
    nothing
end

# function get_signal(path::Vector{CartesianPoint{T}}, charge::T, wp::Interpolations.Extrapolation{T, 3}, S::AbstractCoordinateSystem)::Vector{T} where {T <: SSDFloat}
#     signal::Vector{T} = zeros(T, length(path))
#     add_signal!(signal, path, charge, wp, S)
#     return signal
# end
