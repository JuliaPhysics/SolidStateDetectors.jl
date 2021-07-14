abstract type AbstractWorld{T, ND} end

struct SSDInterval{T <: SSDFloat, L, R, BL, BR} <: IntervalSets.TypedEndpointsInterval{L,R,T}
    left::T
    right::T
end

function get_boundary_types(int::SSDInterval{T, L, R, BL, BR}) where {T <:SSDFloat, L, R, BL, BR}
    return L, R, BL, BR
end

const boundary_condition_mapping = Dict{String, Symbol}(
    "r0" => :r0,
    "inf" => :infinite,
    "infinite" => :infinite,
    "fixed" => :fixed,
    "fix" => :fixed,
    "reflecting" => :reflecting,
    "periodic" => :periodic,
)

struct World{T <: SSDFloat, ND, S} <: AbstractWorld{T, ND} 
    intervals::NTuple{ND, SSDInterval{T}}
end

function World{T, ND, S}(args...) where {T <: SSDFloat, ND, S} 
    return World{T, ND, S}(args)
end


function World(T, dict::Dict, input_units::NamedTuple)::World
    if dict["coordinates"] == "cylindrical"
        CylindricalWorld(T, dict["axes"], input_units)
    elseif dict["coordinates"] == "cartesian"
        CartesianWorld(T, dict["axes"], input_units)
    else
        error("Gridtype must be \"cylindrical\" or \"cartesian\"")
    end
end

function get_periodicity(int::SSDInterval{T, L, R, BL, BR}) where {T <:SSDFloat, L, R, BL, BR}
    return int.right - int.left
end


function get_interval_boundary_types(dict::Dict)
    BL, BR = missing, missing
    if haskey(dict, "boundaries")
        if typeof(dict["boundaries"]) == String
            BL = boundary_condition_mapping[dict["boundaries"]]
            BR = boundary_condition_mapping[dict["boundaries"]]
        else
            BL = boundary_condition_mapping[dict["boundaries"]["left"]]
            BR = boundary_condition_mapping[dict["boundaries"]["right"]]
        end
    end
    if !ismissing(BL)
        if BL == :periodic || BR == :periodic
            if !(BL == :periodic && BR == :periodic) && !((BL == :periodic && BR == :reflecting) || (BL == :reflecting && BR == :periodic))
                throw(ConfigFileError("both or none endings must be \"periodic\" for an interval. Or, \"periodic\" and \"reflecting\" (for periodic behaviour plus mirror symmetry)."))
            end
        end
    end
    return BL, BR
end

function is_periodic_plus_mirror_symmetric(BL::Symbol, BR::Symbol)::Bool
    return (BL == :periodic && BR == :reflecting) || (BL == :reflecting && BR == :periodic)
end

function get_r_SSDInterval(T, dict, input_units::NamedTuple)
    if isa(dict, Dict)
        from::T = 0 
        if "from" in keys(dict) @warn "ConfigFileWarning: \"from\" is not used in r-axis. It is fixed to 0." end
        # to::T = "to" in keys(dict) ? geom_round(ustrip(uconvert(internal_length_unit, T(dict["to"]) * input_units.length))) : throw(ConfigFileError("No \"to\" given for r-axis."))
        to::T = "to" in keys(dict) ? _parse_value(T, dict["to"], input_units.length) : throw(ConfigFileError("No \"to\" given for r-axis."))
        if from < 0 throw(ConfigFileError("left boundary of r-axis cannot be negative.")) end 
        if to < 0 throw(ConfigFileError("right boundary of r-axis cannot be negative.")) end 
        L = :closed
        BL = :r0
        R, BR = :closed, :infinite
        if "boundaries" in keys(dict)
            BR = boundary_condition_mapping[dict["boundaries"]]
        end      
        return SSDInterval{T, L, R, BL, BR}(from, to)
    else
        SSDInterval{T, :closed, :closed, :r0, :infinite}(0, geom_round(ustrip(uconvert(internal_length_unit, T(dict) * input_units.length))))
    end
end
function get_φ_SSDInterval(T, dict::Dict, input_units::NamedTuple)
    if haskey(dict, "phi")
        dp = dict["phi"]
        from::T = "from" in keys(dp) ?  _parse_value(T, dp["from"], input_units.angle) : T(0)
        to::T = "to" in keys(dp) ? _parse_value(T, dp["to"], input_units.angle) : T(2π)
        L = :closed
        R = :open
        BL = :periodic
        BR = :periodic
        cfBL, cfBR = get_interval_boundary_types(dp)
        if !ismissing(cfBL) BL = cfBL; BR = cfBR; end
        if is_periodic_plus_mirror_symmetric(BL, BR) 
            L = :closed; R = :closed;
            BL = :reflecting; BR = :reflecting
        end
        if from == to # 2D
            L = :closed; R = :closed;
            BL = :reflecting; BR = :reflecting
        end
        return SSDInterval{T, L, R, BL, BR}(from, to)
    else
        return SSDInterval{T, :closed, :open, :periodic, :periodic}(T(0), T(2π))
    end
end

function get_cartesian_SSDInterval(T, dict::Dict, input_units::NamedTuple)
    from::T = "from" in keys(dict) ? _parse_value(T, dict["from"], input_units.length) : 0
    to = "to" in keys(dict) ? _parse_value(T, dict["to"], input_units.length) : throw(ConfigFileError("No \"to\" given for z-axis."))
    L = :closed
    R = :closed
    BL, BR = :infinite, :infinite
    cfBL, cfBR = get_interval_boundary_types(dict)
    if !ismissing(cfBL) BL = cfBL; BR = cfBR; end
    if BL == BR == :periodic
        L = :closed; R = :open
    end
    return SSDInterval{T, L, R, BL, BR}(from, to)
end


function CylindricalWorld(T, dict::Dict, input_units::NamedTuple)::World
    r_int = get_r_SSDInterval(T, dict["r"], input_units)
    φ_int = get_φ_SSDInterval(T, dict, input_units)
    z_int = get_cartesian_SSDInterval(T, dict["z"], input_units)
    return World{T, 3, Cylindrical}( r_int, φ_int, z_int )
end


function CartesianWorld(T, dict::Dict, input_units::NamedTuple)::World
    x_int = get_cartesian_SSDInterval(T, dict["x"], input_units)
    y_int = get_cartesian_SSDInterval(T, dict["y"], input_units)
    z_int = get_cartesian_SSDInterval(T, dict["z"], input_units)
    return World{T, 3, Cartesian}( x_int, y_int, z_int )
end

function CartesianWorld(xl::T, xr::T, yl::T, yr::T, zl::T, zr::T)::World where {T <: SSDFloat}
    Δx::T = (xr - xl) / 10
    Δy::T = (yr - yl) / 10
    Δz::T = (zr - zl) / 10
    x_int = SSDInterval{T, :closed, :closed, :infinite, :infinite}(xl - Δx, xr + Δx)
    y_int = SSDInterval{T, :closed, :closed, :infinite, :infinite}(yl - Δy, yr + Δy)
    z_int = SSDInterval{T, :closed, :closed, :infinite, :infinite}(zl - Δz, zr + Δz)
    return World{T, 3, Cartesian}( x_int, y_int, z_int )
end

function CylindricalWorld(r_max::T, zl::T, zr::T)::World where {T <: SSDFloat}
    r_int = SSDInterval{T, :closed, :closed, :r0, :infinite}(T(0), abs(r_max * T(1.1)))
    φ_int = SSDInterval{T, :closed, :open, :periodic, :periodic}(T(0), T(2π))
    Δz::T = zr - zl
    z_int = SSDInterval{T, :closed, :closed, :infinite, :infinite}(zl - (Δz / 10), zr + (Δz / 10))
    return World{T, 3, Cylindrical}( r_int, φ_int, z_int )
end

function World(::Type{Cylindrical}, limits::NTuple{6, T})::World where {T <: SSDFloat}
    return CylindricalWorld(limits[2], limits[5], limits[6])
end
function World(::Type{Cartesian}, limits::NTuple{6, T})::World where {T <: SSDFloat}
    return CartesianWorld(limits...)
end