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




# plotting
_get_plot_r(c::Cone{T, <:Union{T, AbstractInterval{T}}, <:Any, <:Any}) where {T} =
    (_left_radial_interval(c.r),_right_radial_interval(c.r),_left_radial_interval(c.r),_right_radial_interval(c.r))
_get_plot_r(c::Cone{T, <:Tuple, <:Any, <:Any}) where {T} = 
    (_left_radial_interval(c.r[1]),_right_radial_interval(c.r[1]),_left_radial_interval(c.r[2]),_right_radial_interval(c.r[2]))

_get_plot_φ(c::Cone{T, <:Any, Nothing, <:Any}) where {T} = (T(0), T(2π), true)
_get_plot_φ(c::Cone{T, <:Any, <:AbstractInterval, <:Any}) where {T} = (c.φ.left, c.φ.right, false)

_get_plot_z(c::Cone{T}) where {T} = (_left_linear_interval(c.z), _right_linear_interval(c.z))


function get_plot_points(c::Cone{T}) where {T <: AbstractFloat}
    
    plot_points = LineSegments{T}[]
    
    rbotMin::T, rbotMax::T, rtopMin::T, rtopMax::T = _get_plot_r(c)
    φMin::T, φMax::T, φ_is_full_2π::Bool = _get_plot_φ(c)
    zMin::T, zMax::T = _get_plot_z(c)
    
    φrange = range(φMin, φMax, length = 36)
    
    #bottom circle(s)
    for r in [rbotMin, rbotMax]
        if r == 0 continue end
        push!(plot_points, LineSegments{T}([CartesianPoint{T}(r * cos(φ), r * sin(φ), zMin) for φ in φrange]))
    end
    
    #top circle(s)
    for r in [rtopMin, rtopMax]
        if r == 0 continue end
        push!(plot_points, LineSegments{T}([CartesianPoint{T}(r * cos(φ), r * sin(φ), zMax) for φ in φrange]))
    end
    
    #side line(s)
    for φ in (φ_is_full_2π ? T(0) : [φMin, φMax])    
        if rbotMin != 0 || rtopMin != 0
        push!(plot_points, LineSegments{T}([CartesianPoint{T}(rbotMin * cos(φ), rbotMin * sin(φ), zMin), CartesianPoint{T}(rtopMin * cos(φ), rtopMin * sin(φ), zMax)]))      
        end
        push!(plot_points, LineSegments{T}([CartesianPoint{T}(rbotMax * cos(φ), rbotMax * sin(φ), zMin), CartesianPoint{T}(rtopMax * cos(φ), rtopMax * sin(φ), zMax)]))       
    end
    
    #for incomplete φ: lines of cross-sections
    if !φ_is_full_2π
        for φ in [φMin, φMax]
            push!(plot_points, LineSegments{T}([CartesianPoint{T}(rbotMin * cos(φ), rbotMin * sin(φ), zMin), CartesianPoint{T}(rbotMax * cos(φ), rbotMax * sin(φ), zMin)]))
            push!(plot_points, LineSegments{T}([CartesianPoint{T}(rtopMin * cos(φ), rtopMin * sin(φ), zMax), CartesianPoint{T}(rtopMax * cos(φ), rtopMax * sin(φ), zMax)]))
        end
    end
    plot_points
end