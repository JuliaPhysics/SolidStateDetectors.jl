@inline _linear_endpoints(z::Real) = (-z, z)
@inline _linear_endpoints(z::AbstractInterval) = endpoints(z)

@inline _left_radial_interval(r::T) where {T <: Real} = T(0)
@inline _left_radial_interval(r::AbstractInterval) = r.left

@inline _right_radial_interval(r::Real) = r
@inline _right_radial_interval(r::AbstractInterval) = r.right

@inline _extend_number_to_zero_interval(r::T) where {T <: Real} = T(0)..r
@inline _extend_number_to_zero_interval(r::AbstractInterval) = r

@inline get_angular_interval(::Type{T}, α::AbstractInterval{T}) where {T} = α
@inline get_angular_interval(::Type{T}, α::Nothing) where {T} = T(0)..T(2π)

@inline _in_angular_interval_closed(α::Real, α_int::Nothing, tol::Real = 0) = true
@inline _in_angular_interval_closed(α::Real, α_int::AbstractInterval{T}, tol::Real = 0) where {T} = mod(α - (α_int.left-tol), T(2π)) ≤ (α_int.right+tol) - (α_int.left-tol)

@inline function _in_angular_interval_closed(α::T, α_int::Tuple{T,T}; csgtol::T = csg_default_tol(T)) where {T} 
    m = mod(α - α_int[1], T(2π))
    d = (α_int[2] - α_int[1])
    m ≤ d + csgtol 
end
@inline function _in_angular_interval_open(α::T, α_int::Tuple{T,T}; csgtol::T = csg_default_tol(T)) where {T} 
    csgtol < mod(α - α_int[1], T(2π)) < (α_int[2] - α_int[1]) - csgtol 
end


@inline _in_angular_interval_open(α::T, α_int::Nothing) where {T<:Real} = 0 < mod(α, T(2π)) < T(2π)
@inline _in_angular_interval_open(α::Real, α_int::AbstractInterval{T}) where {T} = 0 < mod(α - α_int.left, T(2π)) < α_int.right - α_int.left

_in_angular_interval_union(val::T, α::AbstractInterval, β::AbstractInterval) where {T<:AbstractFloat} =
_in_angular_interval_closed(val, α) || _in_angular_interval_closed(val, β)

function _is_edge_of_angular_interval_union(val::T, α::AbstractInterval, β::AbstractInterval) where {T<:AbstractFloat}
    tol = 10*geom_atol_zero(T)
    _in_angular_interval_union(val+tol, α, β) ⊻ _in_angular_interval_union(val-tol, α, β)
 end

function union_angular_intervals(α::AbstractInterval{T}, β::AbstractInterval{T}) where {T} #if no intersection will return nothing
    tol = 10*geom_atol_zero(T)
    if α == β
        return α
    elseif α.left == α.right && _in_angular_interval_closed(α.left, β)
        return β
    elseif β.left == β.right && _in_angular_interval_closed(β.left, α)
        return α
    elseif α.left ≥ 0 && β.left ≥ 0 && α.right < T(2π) && β.right < T(2π) #less not less or equal to 2pi
        return isempty(α ∩ β) ? nothing : α ∪ β
    else
        edges = 0
        θ1 = T(0)
        θ2 = T(0)
        if _is_edge_of_angular_interval_union(α.left, α, β)
            θ1 = α.left
            edges = edges + 1
        end
        if _is_edge_of_angular_interval_union(α.right, α, β)
            edges == 0 ? θ1 = α.right : θ2 = α.right
            edges = edges + 1
        end
        if _is_edge_of_angular_interval_union(β.left, α, β) && mod(β.left, T(2π)) ≠ mod(α.left, T(2π)) && mod(β.left, T(2π)) ≠ mod(α.right, T(2π))
            if edges == 0
                θ1 = β.left
            elseif edges == 1
                θ2 = β.left
            end
            edges = edges + 1
        end
        if edges < 3 && _is_edge_of_angular_interval_union(β.right, α, β) && mod(β.right, T(2π)) ≠ mod(α.left, T(2π)) && mod(β.right, T(2π)) ≠ mod(α.right, T(2π))
            edges == 1 ? θ2 = β.right : nothing
            edges = edges + 1
        end
        if edges == 2
            θMin, θMax = minmax(mod(θ1, T(2π)), mod(θ2, T(2π)))
            return _in_angular_interval_union(θMin + tol, α, β) ? θMin..θMax : (θMax-T(2π))..θMin
        elseif edges == 0
            return T(0)..T(2π)
        end
    end
end

function is_intersection_an_interval(α::AbstractInterval{T}, β::AbstractInterval{T}) where {T}
    if (mod(α.right - α.left, T(2π)) == 0 || mod(β.right - β.left, T(2π)) == 0)
        return true
    elseif isnothing(union_angular_intervals(α, β))
        return false
    elseif (mod(α.left, T(2π)) == mod(β.right, T(2π)) && !_in_angular_interval_open(α.right, β)) ||
           (mod(α.right, T(2π)) == mod(β.left, T(2π)) && !_in_angular_interval_open(α.left, β))
        return false
    else
        return true
    end
end
