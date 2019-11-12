@inline function get_interpolation(wp::Interpolations.Extrapolation{T, 3}, pt::CartesianPoint{T}, S::Val{:cartesian})::T where {T <: SSDFloat}
    return wp(pt.x, pt.y, pt.z)::T
end
@inline function get_interpolation(wp::Interpolations.Extrapolation{T, 3}, pt::CartesianPoint{T}, S::Val{:cylindrical})::T where {T <: SSDFloat}
    p::CylindricalPoint{T} = CylindricalPoint(pt)
    return wp(p.r, p.φ, p.z)::T
end


function add_signal!(signal::AbstractVector{T}, timestamps::AbstractVector{T}, path::Vector{CartesianPoint{T}}, pathtimestamps::AbstractVector{T}, charge::T, wp::Interpolations.Extrapolation{T, 3}, S::Union{Val{:cylindrical}, Val{:cartesian}})::Nothing where {T <: SSDFloat}
    nsamples::Int = length(pathtimestamps)
    @inbounds for i in eachindex(path)
        signal[i] += get_interpolation(wp, path[i], S)::T * charge
    end
    itp = interpolate( (pathtimestamps,), signal[1:nsamples], Gridded(Linear()))
    signal[1:nsamples] = itp(pathtimestamps)
    if length(signal) > length(path)
        @inbounds signal[length(path)+1:length(signal)] .+= get_interpolation(wp, path[end], S)::T * charge
    end
    nothing
end
function add_signal!(signal::AbstractVector{T}, timestamps::AbstractVector{T}, path::EHDriftPath{T}, charge::T, wp::Interpolations.Extrapolation{T, 3}, S::Union{Val{:cylindrical}, Val{:cartesian}})::Nothing where {T <: SSDFloat}
    add_signal!(signal, timestamps, path.e_path, path.timestamps_e, -charge, wp, S) # electrons induce negative charge
    add_signal!(signal, timestamps, path.h_path, path.timestamps_h,  charge, wp, S)
    nothing
end
function add_signal!(signal::AbstractVector{T}, timestamps::AbstractVector{T}, paths::Vector{EHDriftPath{T}}, charges::Vector{T}, wp::Interpolations.Extrapolation{T, 3}, S::Union{Val{:cylindrical}, Val{:cartesian}})::Nothing where {T <: SSDFloat}
    for ipath in eachindex(paths)
        add_signal!(signal, timestamps, paths[ipath], charges[ipath], wp, S)
    end
    nothing
end

function get_signal(path::Vector{CartesianPoint{T}}, charge::T, wp::Interpolations.Extrapolation{T, 3}, S::Union{Val{:cylindrical}, Val{:cartesian}})::Vector{T} where {T <: SSDFloat}
    signal::Vector{T} = zeros(T, length(path))
    add_signal!(signal, path, charge, wp, S)
    return signal
end
