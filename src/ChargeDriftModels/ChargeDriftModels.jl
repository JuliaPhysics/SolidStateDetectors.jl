abstract type AbstractChargeDriftModel{T <: SSDFloat} end
abstract type AbstractTemperatureModel{T <: SSDFloat} end

"""
    get_electron_drift_field(ef::Array{SVector{3, T},3}, chargedriftmodel::AbstractChargeDriftModel)::Array{SVector{3,T},3} where {T <: SSDFloat}

Applies the charge drift model onto the electric field vectors. The field vectors have to be in cartesian coordinates.
"""
function get_electron_drift_field(ef::Array{SVector{3, T},3}, chargedriftmodel::AbstractChargeDriftModel; 
            use_nthreads::Int = Base.Threads.nthreads())::Array{SVector{3,T},3} where {T <: SSDFloat}
    df = Array{SVector{3,T}, 3}(undef, size(ef))
    @onthreads 1:use_nthreads for i3 in workpart(axes(df, 3), 1:use_nthreads, Base.Threads.threadid())
        for i2 in axes(df, 2)
            for i1 in axes(df, 1)
                @inbounds df[i1, i2, i3] = getVe(ef[i1, i2, i3], chargedriftmodel)
            end
        end
    end
    return df
end

"""
    get_hole_drift_field(ef::Array{SVector{3, T},3}, chargedriftmodel::AbstractChargeDriftModel)::Array{SVector{3,T},3} where {T <: SSDFloat}

Applies the charge drift model onto the hole field vectors. The field vectors have to be in cartesian coordinates.
"""
function get_hole_drift_field(ef::Array{SVector{3,T},3}, chargedriftmodel::AbstractChargeDriftModel; 
            use_nthreads::Int = Base.Threads.nthreads())::Array{SVector{3,T},3} where {T <: SSDFloat}
    df = Array{SVector{3,T}, 3}(undef, size(ef))
    @onthreads 1:use_nthreads for i3 in workpart(axes(df, 3), 1:use_nthreads, Base.Threads.threadid())
        for i2 in axes(df, 2)
            for i1 in axes(df, 1)
                @inbounds df[i1, i2, i3] = getVh(ef[i1, i2, i3], chargedriftmodel)
            end
        end
    end
    return df
end

include("Vacuum/Vacuum.jl")
include("ADL/ADL.jl")