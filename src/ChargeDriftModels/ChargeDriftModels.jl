abstract type AbstractChargeDriftModel{T <: SSDFloat} end
abstract type AbstractTemperatureModel{T <: SSDFloat} end

"""
    get_electron_drift_field(ef::Array{SVector{3, T},3}, chargedriftmodel::AbstractChargeDriftModel)::Array{SVector{3,T},3} where {T <: SSDFloat}

Applies the electron charge drift model onto the electric field vectors and returns the electron drift field.

!!! note 
    The field vectors in `ef` have to be in Cartesian coordinates. 

## Arguments
* `ef::Array{SVector{3, T},3}`: Three-dimensional array with [`ElectricField`](@ref) vectors, e.g. `simulation.electric_field.data`.
* `chargedriftmodel::AbstractChargeDriftModel`: Model that describes the electron drift in the [`Semiconductor`](@ref).

## Keywords 
* `use_nthreads::Int = Base.Threads.nthreads()`: Number of threads that should be used when calculating the drift fields.

## Example 
```julia 
get_electron_drift_field(simulation.electric_field.data, simulation.detector.semiconductor.charge_drift_model)
```

!!! note 
    This method only works if `sim.electric_field` has already been calculated and is not `missing`.
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

Applies the hole charge drift model onto the electric field vectors and returns the hole drift field.

!!! note 
    The field vectors in `ef` have to be in Cartesian coordinates. 

## Arguments
* `ef::Array{SVector{3, T},3}`: Three-dimensional array with [`ElectricField`](@ref) vectors, e.g. `simulation.electric_field.data`.
* `chargedriftmodel::AbstractChargeDriftModel`: Model that describes the electron drift in the [`Semiconductor`](@ref).

## Keywords 
* `use_nthreads::Int = Base.Threads.nthreads()`: Number of threads that should be used when calculating the drift fields.

## Example 
```julia 
get_hole_drift_field(simulation.electric_field.data, simulation.detector.semiconductor.charge_drift_model)
```

!!! note 
    This method only works if `sim.electric_field` has already been calculated and is not `missing`.
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

include("ElectricFieldChargeDriftModel/ElectricFieldChargeDriftModel.jl")
include("ADLChargeDriftModel/ADLChargeDriftModel.jl")