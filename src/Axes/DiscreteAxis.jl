abstract type AbstractAxis{T, BL, BR, I} <: AbstractVector{T} end

"""
    struct DiscreteAxis{T, BL, BR, I} <: AbstractAxis{T, BL, BR, I}
    
Axis with discrete ticks which is used to define a dimension of a [`Grid`](@ref).

## Parametric types 
* `T`: Type of ticks
* `BL`: Boundary condition at the left endpoint.
* `BR`: Boundary condition at the right endpoint.
* `I`: IntervalSets.Interval (closed or open boundaries)

The boundary conditions of a `DiscreteAxis` can be
`BL, BR ∈ {:periodic, :reflecting, :infinite, :r0, :fixed}`.

## Fields
* `interval::I`: Interval that defines the range of the axis.
* `ticks::Vector{T}`: Array of values that correspond to the discrete ticks of the axis.

See also [`Grid`](@ref).
"""
struct DiscreteAxis{T, BL, BR, I} <: AbstractAxis{T, BL, BR, I} 
    interval::I
    ticks::Vector{T}
end
@inline size(dx::DiscreteAxis{T, BL, BR}) where {T, BL, BR} = size(dx.ticks)
@inline IndexStyle(::Type{<:DiscreteAxis}) = IndexLinear()
@inline length(dx::DiscreteAxis{T, BL, BR}) where {T, BL, BR} = length(dx.ticks)
@inline getindex(dx::DiscreteAxis{T, BL, BR}, i::Int) where {T, BL, BR} = dx.ticks[i]
@inline setindex!(dx::DiscreteAxis{T, BL, BR}, v::T, i::Int) where {T, BL, BR} = setindex!(dx.ticks, v, i)
@inline axes(dx::DiscreteAxis{T, BL, BR}) where {T, BL, BR} = axes(dx.ticks)

function DiscreteAxis{T, BL, BR}(int::I, ticks::Vector{T})::DiscreteAxis{T, BL, BR, typeof(int)} where {T, BL, BR, I}
    return DiscreteAxis{T, BL, BR, typeof(int)}( int, ticks )
end

"""
    DiscreteAxis(left_endpoint::T, right_endpoint::T, BL::Symbol, BR::Symbol, L::Symbol, R::Symbol, ticks::AbstractVector{T}) where {T}

Constructor of a `DiscreteAxis`.

## Arguments
* `left_endpoint::T`: Left endpoint of the interval of the `DiscreteAxis`.
* `right_endpoint::T`: Right endpoint of the interval of the `DiscreteAxis`.
* `BL::Symbol`: Boundary condition at the left endpoint.
* `BR::Symbol`: Boundary condition at the right endpoint.
* `L::Symbol`: Boundary type of the left endpoint.
* `R::Symbol`: Boundary type of the right endpoint.
* `ticks::AbstractVector{T}`: Array of values that correspond to the discrete ticks of the axis.

The boundary conditions of a `DiscreteAxis` can be
`BL, BR ∈ {:periodic, :reflecting, :infinite, :r0, :fixed}`.

The boundary types of a `DiscreteAxis` can be `L, R ∈ {:closed, :open}`.

## Examples 
    DiscreteAxis(-2.0, 2.0, :infinite, :infinite, :closed, :closed, collect(-2:0.1:2))
"""
function DiscreteAxis(left_endpoint::T, right_endpoint::T, BL::Symbol, BR::Symbol, L::Symbol, R::Symbol, ticks::AbstractVector{T}) where {T}
    int::Interval{L, R, T} = Interval{L, R, T}( left_endpoint, right_endpoint )
    return DiscreteAxis{T, BL, BR, typeof(int)}( int, ticks )
end


function sizeof(dx::DiscreteAxis{T, BL, BR}) where {T, BL, BR}
    return sizeof(dx.interval) + sizeof(dx.ticks)
end

function print(io::IO, dx::DiscreteAxis{T, BL, BR}) where {T, BL, BR}
    print(io, dx.interval, " - length = ", length(dx))
end
function println(io::IO, dx::DiscreteAxis{T, BL, BR}) where {T, BL, BR}
    println(io, dx.interval)
    println(io, "length = ", length(dx))
end
function show(io::IO, dx::DiscreteAxis{T, BL, BR}) where {T, BL, BR}
    print(io, dx)
end
function show(io::IO, ::MIME"text/plain", dx::DiscreteAxis{T, BL, BR}) where {T, BL, BR}
    show(io, dx)
end

function get_boundary_types(int::Interval{L, R})::Tuple{Symbol, Symbol} where {L, R}
    return L, R
end

function get_boundary_types(ax::DiscreteAxis{T,LB,RB})::NTuple{4, Symbol} where {T, LB, RB}
    return LB, RB, get_boundary_types(ax.interval)...
end

function get_extended_ticks( ax::DiscreteAxis{T, :reflecting, :reflecting} )::Vector{T} where {T}
    ticks_ext::Vector{T} = Array{T}(undef, length(ax.ticks) + 2)
    ticks_ext[2:end-1] = ax.ticks
    set_periodic_boundary_ticks!(ticks_ext, ax.interval)
    return ticks_ext
end
function get_extended_ticks( ax::DiscreteAxis{T, :fixed, :reflecting} )::Vector{T} where {T}
    ticks_ext::Vector{T} = Array{T}(undef, length(ax.ticks) + 2)
    ticks_ext[2:end-1] = ax.ticks
    set_periodic_boundary_ticks!(ticks_ext, ax.interval)
    return ticks_ext
end
function get_extended_ticks( ax::DiscreteAxis{T, :reflecting, :fixed} )::Vector{T} where {T}
    ticks_ext::Vector{T} = Array{T}(undef, length(ax.ticks) + 2)
    ticks_ext[2:end-1] = ax.ticks
    set_periodic_boundary_ticks!(ticks_ext, ax.interval)
    return ticks_ext
end
function get_extended_ticks( ax::DiscreteAxis{T, :periodic, :periodic} )::Vector{T} where {T}
    ticks_ext::Vector{T} = Array{T}(undef, length(ax.ticks) + 2)
    ticks_ext[2:end-1] = ax.ticks
    set_periodic_boundary_ticks!(ticks_ext, ax.interval)
    return ticks_ext
end
function get_extended_ticks( ax::DiscreteAxis{T, :infinite, :infinite} )::Vector{T} where {T}
    ticks_ext::Vector{T} = Array{T}(undef, length(ax.ticks) + 2)
    ticks_ext[2:end-1] = ax.ticks
    Δ::T = 1 * (ticks_ext[end-1] - ticks_ext[2])
    ticks_ext[1] = ticks_ext[2] - Δ
    ticks_ext[end] = ticks_ext[end - 1] + Δ
    return ticks_ext
end
function get_extended_ticks( ax::DiscreteAxis{T, :r0, :infinite} )::Vector{T} where {T}
    ticks_ext::Vector{T} = Array{T}(undef, length(ax.ticks) + 2)
    ticks_ext[2:end-1] = ax.ticks
    ticks_ext[1] = ticks_ext[2] - (ticks_ext[3] - ticks_ext[2])
    Δ::T = 1 * (ticks_ext[end-1] - ticks_ext[2])
    ticks_ext[end] = ticks_ext[end - 1] + Δ
    return ticks_ext
end
function get_extended_ticks( ax::DiscreteAxis{T, :r0, :fixed} )::Vector{T} where {T}
    ticks_ext::Vector{T} = Array{T}(undef, length(ax.ticks) + 2)
    ticks_ext[2:end-1] = ax.ticks
    ticks_ext[1] = ticks_ext[2] - (ticks_ext[3] - ticks_ext[2])
    Δ::T = 1 * (ticks_ext[end-1] - ticks_ext[2])
    ticks_ext[end] = ticks_ext[end - 1] + Δ
    return ticks_ext
end
function get_extended_ticks( ax::DiscreteAxis{T, :r0, :reflecting} )::Vector{T} where {T}
    ticks_ext::Vector{T} = Array{T}(undef, length(ax.ticks) + 2)
    ticks_ext[2:end-1] = ax.ticks
    ticks_ext[1] = ticks_ext[2] - (ticks_ext[3] - ticks_ext[2])
    Δ::T = ticks_ext[end-1] - ticks_ext[end - 2]
    ticks_ext[end] = ticks_ext[end - 1] + Δ
    return ticks_ext
end
function get_extended_ticks( ax::DiscreteAxis{T, :fixed, :fixed} )::Vector{T} where {T}
    # same as get_extended_ticks( ax::DiscreteAxis{T, :reflecting, :reflecting} )::Vector{T} where {T}
    ticks_ext::Vector{T} = Array{T}(undef, length(ax.ticks) + 2)
    ticks_ext[2:end-1] = ax.ticks
    set_periodic_boundary_ticks!(ticks_ext, ax.interval)
    return ticks_ext
end
function get_extended_ticks( ax::DiscreteAxis{T, :infinite, :fixed} )::Vector{T} where {T}
    ticks_ext::Vector{T} = Array{T}(undef, length(ax.ticks) + 2)
    ticks_ext[2:end-1] = ax.ticks
    set_periodic_boundary_ticks!(ticks_ext, ax.interval)
    Δ::T = 1 * (ticks_ext[end-1] - ticks_ext[2])
    ticks_ext[1] = ticks_ext[2] - Δ
    return ticks_ext
end
function get_extended_ticks( ax::DiscreteAxis{T, :infinite, :reflecting} )::Vector{T} where {T}
    ticks_ext::Vector{T} = Array{T}(undef, length(ax.ticks) + 2)
    ticks_ext[2:end-1] = ax.ticks
    set_periodic_boundary_ticks!(ticks_ext, ax.interval)
    Δ::T = 1 * (ticks_ext[end-1] - ticks_ext[2])
    ticks_ext[1] = ticks_ext[2] - Δ
    return ticks_ext
end
function get_extended_ticks( ax::DiscreteAxis{T, :fixed, :infinite} )::Vector{T} where {T}
    ticks_ext::Vector{T} = Array{T}(undef, length(ax.ticks) + 2)
    ticks_ext[2:end-1] = ax.ticks
    set_periodic_boundary_ticks!(ticks_ext, ax.interval)
    Δ::T = 1 * (ticks_ext[end-1] - ticks_ext[2])
    ticks_ext[end] = ticks_ext[end - 1] + Δ
    return ticks_ext
end
function get_extended_ticks( ax::DiscreteAxis{T, :reflecting, :infinite} )::Vector{T} where {T}
    ticks_ext::Vector{T} = Array{T}(undef, length(ax.ticks) + 2)
    ticks_ext[2:end-1] = ax.ticks
    set_periodic_boundary_ticks!(ticks_ext, ax.interval)
    Δ::T = 1 * (ticks_ext[end-1] - ticks_ext[2])
    ticks_ext[end] = ticks_ext[end - 1] + Δ
    return ticks_ext
end

function set_periodic_boundary_ticks!( ticks::Vector{T}, interval::Interval{:closed, :open, T})::Nothing where {T}
    ticks[1] = ticks[2] - (interval.right - ticks[end - 1])
    ticks[end] = interval.right
    nothing
end
function set_periodic_boundary_ticks!( ticks::Vector{T}, interval::Interval{:open, :closed, T})::Nothing where {T}
    ticks[1] = interval.left
    ticks[end] = ticks[end - 1] + (ticks[2] - interval.left)
    nothing
end
function set_periodic_boundary_ticks!( ticks::Vector{T}, interval::Interval{:open, :open, T})::Nothing where {T}
    ticks[1] = interval.left
    ticks[end] = interval.right
    nothing
end

function set_periodic_boundary_ticks!( ticks::Vector{T}, interval::Interval{:closed, :closed, T})::Nothing where {T}
    if length(ticks) == 3 
        ticks[1] = ticks[2] - 2π
        ticks[end] = ticks[2] + 2π # -> Δmidpoint_φ = 2π -> area of circle is 2π * 0.5*r^2   
    else
        ticks[1] = ticks[2] - (ticks[3] - ticks[2])
        ticks[end] = ticks[end - 1] + (ticks[end - 1] - ticks[end - 2])
    end
    nothing
end


function searchsortednearest(a::AbstractVector{T}, x::T)::Int where {T <: Real}
    idx::Int = searchsortedfirst(a, x)
    if (idx == 1) return idx end
    if (idx > length(a)) return length(a) end
    if (a[idx] == x) return idx end
    if (abs(a[idx] - x) < abs(a[idx - 1] - x))
        return idx
    else
        return idx - 1
    end
end
@inline function searchsortednearest(ax::DiscreteAxis{T}, x::T)::Int where {T <: Real}
    return searchsortednearest(ax.ticks, x)
end
function searchsortednearest(ax::DiscreteAxis{T, :periodic, :periodic}, x::T)::Int where {T <: Real}
    if x in ax.interval
        return searchsortednearest(ax.ticks, x)
    else
        period::T = width(ax.interval)
        v::T = x
        while v >= ax.interval.right
            v -= period
        end
        while v < ax.interval.left
            v += period
        end
        return searchsortednearest(ax.ticks, v)
    end    

end

function DiscreteAxis(nt::NamedTuple; unit = Unitful.NoUnits)
    T = typeof(ustrip(nt.knots[1]))
    knots::Vector{T} = convert(Vector{T}, ustrip.(unit, nt.knots))
    lep::T = ustrip(unit, nt.interval.left_boundary.endpoint )
    rep::T = ustrip(unit, nt.interval.right_boundary.endpoint)
    int = Interval{nt.interval.left_boundary.closedopen, nt.interval.right_boundary.closedopen}( lep, rep )
    return DiscreteAxis{T, nt.interval.left_boundary.boundaryhandling, nt.interval.right_boundary.boundaryhandling, typeof(int)}(
        int, knots
    )
end

Base.convert(T::Type{DiscreteAxis}, x::NamedTuple; unit = Unitful.NoUnits) = T(x, unit = unit)

function Base.NamedTuple(ax::DiscreteAxis{T, BL, BR}; unit = Unitful.NoUnits) where {T, BL, BR}
    int::Interval = ax.interval
    int_types::Tuple{Symbol, Symbol} = get_boundary_types(int)
    return (
        knots = ax.ticks * unit, 
        interval = (
            left_boundary = (
                endpoint = int.left * unit,
                closedopen = int_types[1],
                boundaryhandling = BL,
            ),
            right_boundary = (
                endpoint = int.right * unit,
                closedopen = int_types[2],
                boundaryhandling = BR,
            ),
        )
    )
end

Base.convert(T::Type{NamedTuple}, x::DiscreteAxis; unit = Unitful.NoUnits) = T(x, unit = unit)

function merge_closest_ticks!(v::AbstractVector{T}, n::Int = length(v); min_diff::T = T(1e-6)) where {T}
    n == 1 && return n
    Δv = diff(v[1:n])
    Δ_min, Δv_min_indx = findmin(Δv)
    vFirst = v[1]
    vLast  = v[n]
    if Δ_min < min_diff
        v[Δv_min_indx] = (v[Δv_min_indx]+v[Δv_min_indx+1]) / 2
        v[Δv_min_indx+1:end-1] = v[Δv_min_indx+2:end]
        v[1] = vFirst
        n -= 1
        v[n] = vLast
        n
    else
        n
    end
end
function merge_close_ticks(v::AbstractVector{T}; min_diff::T = T(1e-6)) where {T}
    l = length(v)
    l <= 1 && return v
    n = merge_closest_ticks!(v, min_diff = min_diff)
    reduced = n < l
    l = n
    while reduced
        n = merge_closest_ticks!(v, n, min_diff = min_diff)
        reduced = n < l
        l = n
    end
    v[1:n]
end

"""
    merge_second_order_important_ticks(imp::Vector{T}, imp_second_order::Vector{T}; min_diff::T = T(1e-6)) where {T}

Merge all elements of the second vector, `imp_second_order`, into the first vector, `imp`, 
if they are not too close (via `min_diff`) to elements of the first vector.
Returns the merged vector sorted.
"""
function merge_second_order_important_ticks(imp::Vector{T}, imp_second_order::Vector{T}; min_diff::T = T(1e-6)) where {T}
    sorted_imp = sort(imp)
    nearest_inds = map(x -> searchsortednearest(sorted_imp, x), imp_second_order)
    merge_inds = filter(i -> abs(sorted_imp[nearest_inds[i]] - imp_second_order[i]) >= min_diff, eachindex(imp_second_order))
    return sort!(vcat(imp, imp_second_order[merge_inds]))
end

function get_new_ticks_to_equalize_ratios_on_side(t::AbstractVector{T}; max_ratio = T(2)) where {T}
    @assert length(t) == 3 
    Δt = diff(t)
    r = Δt[2] / Δt[1]
    min_ratio = inv(max_ratio)
    new_ticks = T[]
    if r > max_ratio
        x0 = t[2]
        Δnt = max_ratio * Δt[1]
        nt = x0 + Δnt # new tick 
        while t[3] - nt > Δnt
            push!(new_ticks, nt)
            Δnt = max_ratio * Δnt
            nt = nt + Δnt # new tick 
        end
        vcat(t[1:2], new_ticks, t[3:3])
    elseif r < min_ratio
        x0 = t[2]
        Δnt = max_ratio * Δt[2]
        nt = x0 - Δnt # new tick 
        while nt - t[1] > Δnt
            push!(new_ticks, nt)
            Δnt = max_ratio * Δnt
            nt = nt - Δnt # new tick 
        end
        sort!(new_ticks)
        vcat(t[1:1], new_ticks, t[2:3])
    else
        t
    end
end

function initialize_axis_ticks(t::AbstractVector{T}; max_ratio = T(2)) where {T}
    @assert max_ratio >= 1 
    if length(t) <= 2 return t end 
    # First: Left and right side intervals:
    new_ticks_left = get_new_ticks_to_equalize_ratios_on_side(t[1:3]; max_ratio)
    if length(t) == 3 return new_ticks_left end
    if length(t) == 4
        return vcat(new_ticks_left, t[4:4])
    end
    new_ticks_right = get_new_ticks_to_equalize_ratios_on_side(t[end-2:end]; max_ratio)
    if length(t) == 5
        return vcat(new_ticks_left, new_ticks_right[2:end])
    end
    ticks = vcat(new_ticks_left, t[4:end-3], new_ticks_right)
    # Second: Intervals in between
    iL = length(new_ticks_left) - 1 
    iR = iL + 3
    n_sub_ints = length(t) - 5
    for i_sub_int in 1:n_sub_ints
        subticks = ticks[iL:iR]
        Δst = diff(subticks)
        r_l = Δst[2] / Δst[1]
        r_r = Δst[2] / Δst[3]
        r_lr = Δst[1] / Δst[3]
        new_ticks = T[]
        if r_lr < 1 # interval to the left is smaller than interval to the right
            if r_l > max_ratio
                x0 = subticks[2]
                Δnt = max_ratio * Δst[1]
                nt = x0 + Δnt # new tick 
                while subticks[3] - nt > Δnt
                    push!(new_ticks, nt)
                    Δnt = max_ratio * Δnt
                    nt = nt + Δnt # new tick 
                end
            end
        else # interval to the right is smaller than interval to the left
            if r_r > max_ratio
                x0 = subticks[3]
                Δnt = max_ratio * Δst[3]
                nt = x0 - Δnt # new tick 
                while nt - subticks[2] > Δnt
                    push!(new_ticks, nt)
                    Δnt = max_ratio * Δnt
                    nt = nt - Δnt # new tick 
                end
                sort!(new_ticks)
            end
        end
        
        ticks = vcat(ticks[1:iL+1], new_ticks, ticks[iR-1:end])

        iL += 1 + length(new_ticks)
        iR += 1 + length(new_ticks)
    end
    @assert issorted(ticks)
    @assert allunique(ticks)
    ticks
end

function fill_up_ticks(v::AbstractVector{T}, max_diff::T) where {T}
    if length(v) == 1 return v end
    Δv = diff(v)
    add_n_points = Int.(round.(Δv ./ max_diff, RoundUp)) .- 1
    r = Vector{T}(undef, length(v) + sum(add_n_points))
    r[:] .= 0
    i = 0
    for j in eachindex(Δv)
        x0 = v[j]
        n = add_n_points[j]
        Δ = Δv[j] / (n + 1) 
        for l in 0:n
            i += 1
            r[i] = x0 + Δ * l
        end
    end
    r[end] = v[end]
    r
end

function even_tick_axis(ax::DiscreteAxis) 
    if isodd(length(ax)) 
        int = ax.interval
        ticks = ax.ticks
        imax = findmax(diff(ticks))[2]
        push!(ticks, (ticks[imax] + ticks[imax+1]) / 2)
        sort!(ticks)
        typeof(ax)(int, ticks)
    else 
        ax
    end
end

multiplicity(g::DiscreteAxis{T, :infinite, :infinite, I}, ::Type{Cartesian}) where {T, I} = one(T)
multiplicity(g::DiscreteAxis{T, :reflecting, :infinite, I}, ::Type{Cartesian}) where {T, I} = T(2)
multiplicity(g::DiscreteAxis{T, :infinite, :reflecting, I}, ::Type{Cartesian}) where {T, I} = T(2)
multiplicity(g::DiscreteAxis{T, :fixed, :fixed, I}, ::Type{Cartesian}) where {T, I} = one(T)
multiplicity(g::DiscreteAxis{T, :fixed, :reflecting, I}, ::Type{Cartesian}) where {T, I} = T(2)
multiplicity(g::DiscreteAxis{T, :reflecting, :fixed, I}, ::Type{Cartesian}) where {T, I} = T(2)
multiplicity(g::DiscreteAxis{T, :fixed, :infinite, I}, ::Type{Cartesian}) where {T, I} = one(T)
multiplicity(g::DiscreteAxis{T, :infinite, :fixed, I}, ::Type{Cartesian}) where {T, I} = one(T)

function multiplicity(g::DiscreteAxis{T, :reflecting, :reflecting, I}, ::Type{Cartesian}) where {T, I} 
    @warn "Multiplicity of Cartesian axis, $(g) (:reflecting, :reflecting), would be infinite. It is set to 1 here." 
    one(T)
end
