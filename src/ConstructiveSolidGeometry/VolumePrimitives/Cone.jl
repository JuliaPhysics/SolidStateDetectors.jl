struct Cone{T,TR,TP,TZ} <: AbstractVolumePrimitive{T}
    r::TR #if not a Tuple, then Cone is a Tube
    φ::TP
    z::TZ
    #if r is a Tuple, the first entry refers to the r-interval at the bottom, the second one to the r-interval at the top
    function Cone( ::Type{T},
                   r::Union{T, <:AbstractInterval{T}, Tuple{T,T}, Tuple{I,I}},
                   φ::Union{Nothing, <:AbstractInterval{T}},
                   z::Union{T, <:AbstractInterval{T}}) where {T, I<:AbstractInterval{T}}
        new{T,typeof(r),typeof(φ),typeof(z)}(r, φ, z)
    end
end

#Constructors
function Cone(;rbotMin = 0, rbotMax = 1, rtopMin = 0, rtopMax = 1, φMin = 0, φMax = 2π, zMin = -1/2, zMax = 1/2)
    T = float(promote_type(typeof.((rtopMin, rtopMax, rbotMin, rbotMax, φMin, φMax, zMin, zMax))...))
    rMin_is_equal::Bool = rbotMin == rtopMin
    rMax_is_equal::Bool = rbotMax == rtopMax
    rMin_is_zero::Bool = rMin_is_equal && rbotMin == 0
    r = if rMax_is_equal
            if rMin_is_zero # Tube with rMin = 0
                T(rbotMax) 
            elseif rMin_is_equal # Tube
                T(rbotMin)..T(rbotMax)
            else # Cone
                (T(rbotMin)..T(rbotMax), T(rtopMin)..T(rtopMax))
            end
        elseif rMin_is_zero #Cone with rMin = 0
            (T(rbotMax), T(rtopMax))
        else # Cone
            (T(rbotMin)..T(rbotMax), T(rtopMin)..T(rtopMax))
        end
    φ = mod(T(φMax) - T(φMin), T(2π)) == 0 ? nothing : T(φMin)..T(φMax)
    z = zMax == -zMin ? T(zMax) : T(zMin)..T(zMax)
    Cone( T, r, φ, z)
end
Cone(rbotMin, rbotMax, rtopMin, rtopMax, φMin, φMax, zMin, zMax) = Cone(;rbotMin, rbotMax, rtopMin, rtopMax, φMin, φMax, zMin, zMax)

function Cone(rbot::R1, rtop::R2, height::H) where {R1<:Real, R2<:Real, H<:Real}
    T = float(promote_type(R1, R2, H))
    Cone( T, (T(rbot), T(rtop)), nothing, T(height)/2)
end


#Constructors for Tubes
Tube(;rMin = 0, rMax = 1, φMin = 0, φMax = 2π, zMin = -1/2, zMax = 1/2) = Cone(rMin, rMax, rMin, rMax, φMin, φMax, zMin, zMax)
Tube(rMin, rMax, φMin, φMax, zMin, zMax) = Tube(;rMin, rMax, φMin, φMax, zMin, zMax) 

function Tube(r::R, height::H) where {R<:Real, H<:Real}
    T = float(promote_type(R,H))
    Cone(T, T(r), nothing, T(height)/2)
end

function Tube(rMin::R1, rMax::R2, height::H) where {R1<:Real, R2<:Real, H<:Real}
    T = float(promote_type(R1,R2,H))
    Cone(T, rMin == 0 ? T(rMax) : T(rMin)..T(rMax), nothing, T(height)/2)
end


# for Tubes
get_r_at_z(c::Cone{T, <:Union{T, AbstractInterval{T}}, <:Any, <:Any}, z::Real) where {T} = c.r 

# for Cones
get_r_at_z(c::Cone{T, Tuple{T,T}, <:Any, <:Any}, z::Real) where {T} = _get_r_at_z(c.r[1], c.r[2], c.z, z)

function get_r_at_z(c::Cone{T, Tuple{I,I}, <:Any, <:Any}, z::Real) where {T, I<:AbstractInterval{T}}
    r1::T = _get_r_at_z(c.r[1].left, c.r[2].left, c.z, z)
    r2::T = _get_r_at_z(c.r[1].right, c.r[2].right, c.z, z)
    r1..r2
end

function _get_r_at_z(rbot::TR, rtop::TR, cz::TZ, z::Real)::TR where {TR<:Real, TZ} 
    (rtop - rbot) * (z - _left_linear_interval(cz)) / _width_linear_interval(cz) + rbot
end


in(p::AbstractCoordinatePoint, c::Cone{<:Any, <:Any, Nothing, <:Any}) =
    _in_z(p, c.z) && _in_cyl_r(p, get_r_at_z(c, p.z))

in(p::AbstractCoordinatePoint, c::Cone{<:Any, <:Any, <:AbstractInterval, <:Any}) =
    _in_z(p, c.z) && _in_φ(p, c.φ) && _in_cyl_r(p, get_r_at_z(c, p.z))
