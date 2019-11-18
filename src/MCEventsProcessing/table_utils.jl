
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
    if all(isequal(row_lengths[1]), row_lengths)
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


function cluster_detector_hits(
    detno::AbstractVector{<:Integer},
    edep::AbstractVector{<:RealQuantity},
    pos::AbstractVector{<:StaticVector{3,PT}},
    cluster_radius::RealQuantity
) where {TT<:RealQuantity,PT <: RealQuantity}
    Table = TypedTables.Table
    unsorted = Table(detno = detno, edep = edep, pos = pos)
    sorting_idxs = sortperm(unsorted.detno)
    sorted = unsorted[sorting_idxs]
    grouped = Table(consgroupedview(sorted.detno, Tables.columns(sorted)))

    r_detno = similar(detno, 0)
    r_edep = similar(edep, 0)
    r_pos = similar(pos, 0)

    posunit = unit(PT)
    ustripped_cradius = ustrip(uconvert(posunit, cluster_radius))
    
    for d_hits_nt in grouped
        d_hits = Table(d_hits_nt)
        d_detno = first(d_hits.detno)
        @assert all(isequal(d_detno), d_hits.detno)
        if length(d_hits) > 3
            clusters = Clustering.dbscan(ustrip.(flatview(d_hits.pos)), ustripped_cradius, leafsize = 20, min_neighbors = 1, min_cluster_size = 1)
            for c in clusters
                idxs = vcat(c.boundary_indices, c.core_indices)
                @assert length(idxs) == c.size
                c_hits = view(d_hits, idxs)
                global c_hits = c_hits
                
                push!(r_detno, d_detno)
                push!(r_edep, sum(c_hits.edep))
                esum = ustrip(r_edep[end])
                inv_e_sum = inv(esum)
                weights = Weights(ustrip.(c_hits.edep), esum) .* inv_e_sum
                push!(r_pos, sum(c_hits.pos .* weights))
            end
        else
            append!(r_detno, d_hits.detno)
            append!(r_edep, d_hits.edep)
            append!(r_pos, d_hits.pos)
        end
    end

    (detno = r_detno, edep = r_edep, pos = r_pos)
end


function cluster_detector_hits(table::TypedTables.Table, cluster_radius::RealQuantity)
    @assert :pos in TypedTables.columnnames(table) "Table has no column `pos`"
    @assert :edep in TypedTables.columnnames(table) "Table has no column `edep`"
    @assert :detno in TypedTables.columnnames(table) "Table has no column `detno`"
    clustered_nt = map(
        evt -> cluster_detector_hits(evt.detno, evt.edep, evt.pos, cluster_radius),
        table
    )
    TypedTables.Table(merge(
        TypedTables.columns(table),
        map(
            VectorOfVectors,
            TypedTables.columns(clustered_nt)
        )
    ))
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
