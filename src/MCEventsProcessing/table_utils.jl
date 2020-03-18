
@inline remove_column(table::TypedTables.Table, colname::Symbol) =
    TypedTables.Table(Base.structdiff(TypedTables.columns(table), NamedTuple{(colname,),Tuple{Int}}(0)))

@inline add_column(table::TypedTables.Table, nt::NamedTuple) = 
    TypedTables.Table(merge(nt, TypedTables.columns(table)))
@inline add_column(table::TypedTables.Table, newname::Symbol, col::CT) where {CT} =
    add_column(table, NamedTuple{(newname,), Tuple{CT}}((col,)) )

@inline _rename_fields_impl(nt::NamedTuple{names}, from::Val{a}, to::Val{b}) where {names, a, b} =
    NamedTuple{map(x -> x == a ? b : x, names)}(values(nt))

@inline rename_fields(nt::NamedTuple, pattern::Pair{Symbol,Symbol}) =
    _rename_fields_impl(nt, Val(pattern.first), Val(pattern.second))

@inline rename_cols(table::Any, pattern::Pair{Symbol,Symbol}) =
    TypedTables.Tables.materializer(table)(rename_fields(TypedTables.Tables.columns(table), pattern))


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

function blow_up_to_commen_length(table::TypedTables.Table)
    n = common_length.(values(TypedTables.columns(table))...)
    return TypedTables.Table(map(c -> nested_col(c, n), TypedTables.columns(table)))
end

function split_table_by_each_charge_deposition(table::TypedTables.Table)
    row_lengths = length.(values(table[1]))
    all_vectors = all(t -> t <: AbstractVector, typeof.(values(table[1])))
    if all(isequal(row_lengths[1]), row_lengths) && all_vectors
        TypedTables.Table( map(flatview, TypedTables.columns(table)) )
    else
        split_table_by_each_charge_deposition(blow_up_to_commen_length(table))
    end
end

function sort_by_time(table::TypedTables.Table)  
    @assert length(table.thit[1]) == 1 "Split the table first: split_table_by_each_charge_deposition"
    sorted_inds = sortperm(table.thit) 
    TypedTables.Table(map( c -> c[sorted_inds], TypedTables.columns(table))) 
end

function cluster_by_time(table::TypedTables.Table, Δt = 1u"ns")
    @assert length(table.thit[1]) == 1 "Split the table first: split_table_by_each_charge_deposition"
    if Δt == 0u"s" return table end
    times_sorted = !issorted(table.thit) ? sort_by_time(table).thit : table.thit
    n_charge_depos = length(times_sorted)
    t0 = times_sorted[1]
    t_last = t0
    inds = Int[]
    for i in 2:n_charge_depos
        t = times_sorted[i]
        dt = t - t_last
        if dt > Δt
            push!(inds, i - 1)
            t_last = t
        end
    end
    if length(inds) > 0
        ranges = vcat( [1:inds[1]], broadcast(:, inds[1:end-1].+1, inds[2:end]), 
                        n_charge_depos > inds[end] ? [inds[end]+1:n_charge_depos] : [] )
        TypedTables.Table(map(c -> VectorOfVectors(map(r -> c[r], ranges)), TypedTables.columns(table)))
    else
        TypedTables.Table(map(c -> VectorOfVectors(map(r -> c[r], [1:n_charge_depos])), TypedTables.columns(table)))
    end
end

function add_physics_event_numbers(table::TypedTables.Table)
    n_events = length(TypedTables.columns(table)[1])
    add_column(table, :phyevtno, flatview(collect(Int32, 1:n_events)))
end

function groub_by_evt_number(table::TypedTables.Table, colname::Symbol = :mcevtno)
    sorting_idxs = sortperm(getproperty(table, colname))
    sorted = table[sorting_idxs]
    return TypedTables.Table(consgroupedview((getproperty(table, colname), Tables.columns(sorted))))
end


function translate_event_positions(table::TypedTables.Table, tv::AbstractVector{<:Unitful.Quantity{<:Real, Unitful.𝐋}})
    @assert length(tv) == 3 "translation vector must be of length 3"
    TP = eltype(eltype(eltype(table.pos)))
    u = unit(TP)
    @assert all(e -> dimension(e) == Unitful.𝐋, tv) "All elements of translation vector muss be of dimension Unitful.𝐋"
    trans = uconvert.(u, tv)
    newpos = (pos = VectorOfVectors(map(poss -> map(pos -> pos + trans, poss), table.pos)),)
    TypedTables.Table(merge( TypedTables.columns(table), newpos ))
end

function rotate_event_positions(table, rot::Rotations.Rotation{3})
    newpos = (pos = VectorOfVectors(map(poss -> map(pos -> rot * pos, poss), table.pos)),)
    TypedTables.Table(merge( TypedTables.columns(table), newpos ))
end

function chunked_ranges(n, chunksize)
    if n < chunksize return [1:n] end
    ranges = map( s -> s:s+chunksize-1, 1:chunksize:n)
    if last(ranges[end]) > n ranges[end] = first(ranges[end]):n end
    return ranges
end
