
@inline add_column(table::TypedTables.Table, nt::NamedTuple) =
    TypedTables.Table(merge(nt, TypedTables.columns(table)))
@inline add_column(table::TypedTables.Table, newname::Symbol, col::CT) where {CT} =
    add_column(table, NamedTuple{(newname,), Tuple{CT}}((col,)) )


function common_length(x, y)
    a = length(x); b = length(y)
    if a == b || b == 1
        a
    elseif a == 1
        b
    else
        throw(ArgumentError("Lengths are not compatible"))
    end
end
common_length(x, y, z...) = common_length(common_length(x, y), z...)

nested_col(x::AbstractVector, n::AbstractVector{<:Integer}) = VectorOfVectors(Fill.(x, n))

function nested_col(x::AbstractVector{<:AbstractVector}, n::AbstractVector{<:Integer})
    @assert length.(x) == n
    x
end

function blow_up_to_common_length(table::TypedTables.Table)
    n = common_length.(values(TypedTables.columns(table))...)
    return TypedTables.Table(map(c -> nested_col(c, n), TypedTables.columns(table)))
end

function split_table_by_each_charge_deposition(table::TypedTables.Table)
    row_lengths = length.(values(table[1]))
    all_vectors = all(t -> t <: AbstractVector, typeof.(values(table[1])))
    if all(isequal(row_lengths[1]), row_lengths) && all_vectors
        TypedTables.Table( map(flatview, TypedTables.columns(table)) )
    else
        split_table_by_each_charge_deposition(blow_up_to_common_length(table))
    end
end

function split_by_time(table::TypedTables.Table, Δt = 1u"ns")
    
    hasproperty(table, :thit) || throw(ArgumentError("Expected detector hit events table to have column named `thit`"))
    elem_ptr = deepcopy(table.thit.elem_ptr)
    t_hit = flatview(table.thit)
    idx = Vector{Int}(undef, length(t_hit))
    for e in 1:size(elem_ptr,1)-1
        t_row = t_hit[elem_ptr[e]:elem_ptr[e+1]-1]
        s = sortperm(t_row)
        t_last = t_row[first(s)]
        for (i,idx) in enumerate(s)
            t = t_row[idx]
            dt = t - t_last
            if _parse_value(Float64, dt, u"s") > _parse_value(Float64, Δt, u"s")
                push!(elem_ptr, i + elem_ptr[e] - 1)
                t_last = t
            end
        end
        idx[elem_ptr[e]:elem_ptr[e+1]-1] .= s .+ (elem_ptr[e] - 1)
    end
    sort!(elem_ptr)
    
    return TypedTables.Table(
        map(c -> all(length.(c) .== 1) ?
            # columns that have a single entry: clone them
            c[map(p -> findlast(p .>= table.thit.elem_ptr), elem_ptr[1:end-1])] :
            # columns that have multiple entries: split accordingly
            VectorOfVectors(flatview(c)[idx], elem_ptr), 
        TypedTables.columns(table))
    )
end

@inline _translate_event_position(p::P, v::AbstractVector{<:Union{<:Real, <:LengthQuantity}}) where {P <: Union{<:CartesianPoint, AbstractVector{<:Real}}} = P(p + to_internal_units(v))
@inline _translate_event_position(p::P, v::AbstractVector{<:Real}) where {P <: AbstractVector{<:LengthQuantity}} = P(p + typeof(ustrip.(p))(v) * internal_length_unit)
@inline _translate_event_position(p::P, v::AbstractVector{<:LengthQuantity}) where {P <: AbstractVector{<:LengthQuantity}} = P(p + P(v))
function translate_event_positions(table::TypedTables.Table, tv::AbstractVector{<:Union{<:Real, <:LengthQuantity}})
    hasproperty(table, :pos) || throw(ArgumentError("Expected detector hit events table to have column named `pos`"))
    length(tv) == 3 || throw(ArgumentError("Translation vector must be of length 3"))
    eltype(eltype(eltype(table.pos))) <: Union{<:Real, <:LengthQuantity} || throw(ArgumentError("Expected table positions to be unitless or to have units of length"))
    newpos = (pos = VectorOfVectors(map(poss -> map(pos -> _translate_event_position(pos, tv), poss), table.pos)),)
    TypedTables.Table(merge( TypedTables.columns(table), newpos ))
end

@inline _rotate_event_position(p::P, rot::Rotations.Rotation{3}) where {P <: CartesianPoint} = P(cartesian_zero + rot * (p - cartesian_zero))
@inline _rotate_event_position(p::P, rot::Rotations.Rotation{3}) where {P <: AbstractVector} = P(rot * p)
function rotate_event_positions(table, rot::Rotations.Rotation{3}) 
    hasproperty(table, :pos) || throw(ArgumentError("Expected detector hit events table to have column named `pos`"))
    newpos = (pos = VectorOfVectors(map(poss -> map(pos -> _rotate_event_position(pos, rot), poss), table.pos)),)
    TypedTables.Table(merge( TypedTables.columns(table), newpos ))
end

function chunked_ranges(n, chunksize)
    if n < chunksize return [1:n] end
    ranges = map( s -> s:s+chunksize-1, 1:chunksize:n)
    if last(ranges[end]) > n ranges[end] = first(ranges[end]):n end
    return ranges
end
