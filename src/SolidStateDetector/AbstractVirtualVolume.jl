abstract type AbstractVirtualVolume{T} end

in(p::AbstractCoordinatePoint{T, 3}, avv::AbstractVirtualVolume{T}) where {T <: SSDFloat} = in(p, avv.geometry)
