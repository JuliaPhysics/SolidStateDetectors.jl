abstract type AbstractChargeCloud end

struct PointCharge{T <: SSDFloat, S} <: AbstractChargeCloud 
    charge::T
    pos::AbstractCoordinatePoint{T, 3, S}
end

struct NBodyChargeCloud{T <: SSDFloat, S} <: AbstractChargeCloud 
    # To be done...
end