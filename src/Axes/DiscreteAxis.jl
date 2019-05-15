abstract type AbstractAxis{T, BL, BR} <: AbstractVector{T} end
abstract type AbstractDiscreteAxis{T, BL, BR <: Number} <: AbstractAxis{T, BL, BR} end

"""
    DiscreteAxis{T, BL, BR} <: AbstractAxis{T, BL, BR}

* T: Type of ticks
* BL, BR ∈ {:periodic, :reflecting, :infinite, :r0, :fixed} 
* BL: left boundary condition
* BR: right boundary condition
"""
struct DiscreteAxis{T, BL, BR} <: AbstractAxis{T, BL, BR} 
    interval::Interval{L, R, T} where {L, R}
    ticks::Vector{T}
end
@inline size(dx::DiscreteAxis{T, BL, BR}) where {T, BL, BR} = size(dx.ticks)
@inline IndexStyle(::Type{<:DiscreteAxis}) = IndexLinear()
@inline length(dx::DiscreteAxis{T, BL, BR}) where {T, BL, BR} = length(dx.ticks)
@inline getindex(dx::DiscreteAxis{T, BL, BR}, i::Int) where {T, BL, BR} = dx.ticks[i]
@inline setindex!(dx::DiscreteAxis{T, BL, BR}, v::T, i::Int) where {T, BL, BR} = setindex!(dx.ticks, v, i)
@inline axes(dx::DiscreteAxis{T, BL, BR}) where {T, BL, BR} = axes(dx.ticks)

"""
    DiscreteAxis(left_endpoint::T, right_endpoint::T, BL::Symbol, BR::Symbol, L::Symbol, R::Symbol, ticks::AbstractVector{T}) where {T}

* T: Type of ticks
* BL, BR ∈ {:periodic, :reflecting, :infinite, :r0, :fixed} 
* L, R {:closed, :open} 
* ticks: Ticks of the axis
"""
function DiscreteAxis(left_endpoint::T, right_endpoint::T, BL::Symbol, BR::Symbol, L::Symbol, R::Symbol, ticks::AbstractVector{T}) where {T}
    int::Interval{L, R, T} = Interval{L, R, T}( left_endpoint, right_endpoint )
    return DiscreteAxis{T, BL, BR}( int, ticks )
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
    println(io, dx)
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


function uniq(v::Vector{T})::Vector{T} where {T <: Real}
    v1::Vector{T} = Vector{T}()
    if length(v) > 0
        laste::T = v[1]
        push!(v1, laste)
        for e in v
            if e != laste
                laste = e
                push!(v1, laste)
            end
        end
    end
    return v1
end

function merge_axis_ticks_with_important_ticks(ax::DiscreteAxis{T}, impticks::Vector{T}; atol::Real = 0.0001 )::Vector{T} where {T}
    v::Vector{T} = T[]
    for r in impticks if in(r, ax.interval) push!(v, r) end end
    for r in ax push!(v, r) end
    sort!(v)
    v = uniq(v)
    delete_idcs::Vector{Int} = Int[]
    for i in 1:(length(v) - 1)
        if (v[i + 1] - v[i]) < atol
            if !in(v[i], impticks) push!(delete_idcs, i) end
            if !in(v[i + 1], impticks) push!(delete_idcs, i + 1) end
        end
    end
    delete_idcs = sort(uniq(delete_idcs))
    deleteat!(v, delete_idcs) 
    for impv in impticks
        if !in(impv, v) && in(impv, ax.interval)
            error("Important ticks were removed.")
        end
    end
    return v
end


function range(interval::Interval{:closed, :closed, T}; step::Union{Missing, T} = missing, length::Union{Missing, Int} = missing) where {T}
    stop::T = interval.right
    if ismissing(step) && ismissing(length)
        range(interval.left, stop = stop)
    elseif ismissing(step)
        range(interval.left, stop = stop, length=length)
    elseif ismissing(length)
        range(interval.left, stop = stop, step=step)
    else
        error(KeyError, ": Both keywords `step` and `length` were given. But only one is allowed.")
    end
end

function range(interval::Interval{:closed, :open, T}; step::Union{Missing, T} = missing, length::Union{Missing, Int} = missing) where {T}
    if ismissing(step) && ismissing(length)
        length::Int = 2
        stop::T = (interval.right + interval.left) / 2
        range(interval.left, stop = stop, length=2)
    elseif ismissing(step)
        stop = interval.right - ( interval.right - interval.left ) / length
        range(interval.left, stop = stop, length=length)
    elseif ismissing(length)
        # stop = interval.right - interval.right % step
        stop = geom_round(interval.right - step)
        range(interval.left, stop = stop, step=step)
    else
        error(KeyError, ": Both keywords `step` and `length` were given. But only one is allowed.")
    end
end


function range(interval::Interval{:open, :closed, T}; step::Union{Missing, T} = missing, length::Union{Missing, Int} = missing) where {T}
    stop::T = interval.right
    if ismissing(step) && ismissing(length)
        step::T = (stop - interval.left) / 2 
        range(interval.left + step, stop = stop, length=2)
    elseif ismissing(step)
        step = (stop - interval.left) / length
        range(interval.left + step, stop = stop, length=length)
    elseif ismissing(length)
        range(interval.left + step, stop = stop, step=step)
    else
        error(KeyError, ": Both keywords `step` and `length` were given. But only one is allowed.")
    end
end

function range(interval::Interval{:open, :open, T}; step::Union{Missing, T} = missing, length::Union{Missing, Int} = missing) where {T}
    if ismissing(step) && ismissing(length)
        step::T = ( interval.right - interval.left ) / 3
        range(interval.left + step, stop = interval.right - step, length=2)
    elseif ismissing(step)
        step = ( interval.right - interval.left ) / (length + 1)
        range(interval.left + step, stop = interval.right - step, length=length)
    elseif ismissing(length)
        tmp::T = interval.right % step
        if tmp == 0 tmp = step end
        stop = interval.right - tmp
        range(interval.left + step, stop = stop, step=step)
    else
        error(KeyError, ": Both keywords `step` and `length` were given. But only one is allowed.")
    end
end

function DiscreteAxis{BL, BR}(interval::Interval{L, R, T}; step::Union{Missing, T} = missing, length::Union{Missing, Int} = missing)::DiscreteAxis{T, BL, BR} where {L, R, T, BL, BR}
    ticks::Vector{T} = collect(range(interval, step=step, length=length))
    if T == Float32 || T == Float64
        ticks = round.(ticks, sigdigits = geom_sigdigits(T))
        for iv in eachindex(ticks)
            if isapprox(ticks[iv], 0, atol = geom_atol_zero(T)) 
                ticks[iv] = zero(T)
            end
        end
    end
    DiscreteAxis{T, BL, BR}(interval, ticks)
end

function midpoints(a::Vector{T})::Vector{T} where {T}
    @inbounds r::Vector{T} = a[1:end-1]
    @simd for i in eachindex(r)
        @inbounds r[i] += 0.5 * (a[i + 1] - a[i])
    end
    return r
end

function get_extended_ticks( ax::DiscreteAxis{T, :reflecting, :reflecting} )::Vector{T} where {T}
    ticks_ext::Vector{T} = Array{T}(undef, length(ax.ticks) + 2)
    ticks_ext[2:end-1] = ax.ticks
    set_periodic_bondary_ticks!(ticks_ext, ax.interval)
    return ticks_ext
end
function get_extended_ticks( ax::DiscreteAxis{T, :fixed, :reflecting} )::Vector{T} where {T}
    ticks_ext::Vector{T} = Array{T}(undef, length(ax.ticks) + 2)
    ticks_ext[2:end-1] = ax.ticks
    set_periodic_bondary_ticks!(ticks_ext, ax.interval)
    return ticks_ext
end
function get_extended_ticks( ax::DiscreteAxis{T, :reflecting, :fixed} )::Vector{T} where {T}
    ticks_ext::Vector{T} = Array{T}(undef, length(ax.ticks) + 2)
    ticks_ext[2:end-1] = ax.ticks
    set_periodic_bondary_ticks!(ticks_ext, ax.interval)
    return ticks_ext
end
function get_extended_ticks( ax::DiscreteAxis{T, :periodic, :periodic} )::Vector{T} where {T}
    ticks_ext::Vector{T} = Array{T}(undef, length(ax.ticks) + 2)
    ticks_ext[2:end-1] = ax.ticks
    set_periodic_bondary_ticks!(ticks_ext, ax.interval)
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
function get_extended_ticks( ax::DiscreteAxis{T, :fixed, :fixed} )::Vector{T} where {T}
    # same as get_extended_ticks( ax::DiscreteAxis{T, :reflecting, :reflecting} )::Vector{T} where {T}
    ticks_ext::Vector{T} = Array{T}(undef, length(ax.ticks) + 2)
    ticks_ext[2:end-1] = ax.ticks
    set_periodic_bondary_ticks!(ticks_ext, ax.interval)
    return ticks_ext
end
function get_extended_ticks( ax::DiscreteAxis{T, :infinite, :fixed} )::Vector{T} where {T}
    ticks_ext::Vector{T} = Array{T}(undef, length(ax.ticks) + 2)
    ticks_ext[2:end-1] = ax.ticks
    set_periodic_bondary_ticks!(ticks_ext, ax.interval)
    Δ::T = 1 * (ticks_ext[end-1] - ticks_ext[2])
    ticks_ext[1] = ticks_ext[2] - Δ
    return ticks_ext
end
function get_extended_ticks( ax::DiscreteAxis{T, :infinite, :reflecting} )::Vector{T} where {T}
    ticks_ext::Vector{T} = Array{T}(undef, length(ax.ticks) + 2)
    ticks_ext[2:end-1] = ax.ticks
    set_periodic_bondary_ticks!(ticks_ext, ax.interval)
    Δ::T = 1 * (ticks_ext[end-1] - ticks_ext[2])
    ticks_ext[1] = ticks_ext[2] - Δ
    return ticks_ext
end
function get_extended_ticks( ax::DiscreteAxis{T, :fixed, :infinite} )::Vector{T} where {T}
    ticks_ext::Vector{T} = Array{T}(undef, length(ax.ticks) + 2)
    ticks_ext[2:end-1] = ax.ticks
    set_periodic_bondary_ticks!(ticks_ext, ax.interval)
    Δ::T = 1 * (ticks_ext[end-1] - ticks_ext[2])
    ticks_ext[end] = ticks_ext[end - 1] + Δ
    return ticks_ext
end
function get_extended_ticks( ax::DiscreteAxis{T, :reflecting, :infinite} )::Vector{T} where {T}
    ticks_ext::Vector{T} = Array{T}(undef, length(ax.ticks) + 2)
    ticks_ext[2:end-1] = ax.ticks
    set_periodic_bondary_ticks!(ticks_ext, ax.interval)
    Δ::T = 1 * (ticks_ext[end-1] - ticks_ext[2])
    ticks_ext[end] = ticks_ext[end - 1] + Δ
    return ticks_ext
end

function set_periodic_bondary_ticks!( ticks::Vector{T}, interval::Interval{:closed, :open, T})::Nothing where {T}
    ticks[1] = ticks[2] - (interval.right - ticks[end - 1])
    ticks[end] = interval.right
    nothing
end
function set_periodic_bondary_ticks!( ticks::Vector{T}, interval::Interval{:open, :closed, T})::Nothing where {T}
    ticks[1] = interval.left
    ticks[end] = ticks[end - 1] + (ticks[2] - interval.left)
    nothing
end
function set_periodic_bondary_ticks!( ticks::Vector{T}, interval::Interval{:open, :open, T})::Nothing where {T}
    ticks[1] = interval.left
    ticks[end] = interval.right
    nothing
end

function set_periodic_bondary_ticks!( ticks::Vector{T}, interval::Interval{:closed, :closed, T})::Nothing where {T, ispolaraxis}
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
        period::T = ax.interval.right - ax.interval.left
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

function DiscreteAxis(nt::NamedTuple; unit = u"m/m")
    T = typeof(ustrip(nt.knots[1]))
    knots::Vector{T} = convert(Vector{T}, ustrip.(uconvert.(unit, nt.knots)))
    lep::T = ustrip(uconvert.(unit, nt.interval.left_boundary.endpoint ))
    rep::T = ustrip(uconvert.(unit, nt.interval.right_boundary.endpoint))
    return DiscreteAxis{T, nt.interval.left_boundary.boundaryhandling, nt.interval.right_boundary.boundaryhandling}(
        Interval{nt.interval.left_boundary.closedopen, nt.interval.right_boundary.closedopen}( lep, rep ),
        knots
    )
end

Base.convert(T::Type{DiscreteAxis}, x::NamedTuple; unit = u"m/m") = T(x)

function NamedTuple(ax::DiscreteAxis{T, BL, BR}; unit = u"m/m") where {T, BL, BR}
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

Base.convert(T::Type{NamedTuple}, x::DiscreteAxis; unit = u"m/m") = T(x)

