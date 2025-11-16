# # This file is a part of SolidStateDetectors.jl, licensed under the MIT License (MIT).


@inline function finalize_group!(time_groups, current_group, new_index)
    push!(time_groups, current_group)
    return [new_index]
end

function cluster_detector_hits(
    detno::AbstractVector{<:Integer},
    edep::AbstractVector{TT},
    pos::AbstractVector{CartesianPoint{PT}},
    thit::AbstractVector{TTT},
    cluster_radius::RealQuantity,
    cluster_time::RealQuantity
) where {TT<:RealQuantity, PT <: RealQuantity, TTT <: RealQuantity}

    pos = [SVector(p.x, p.y, p.z) for p in pos]   # converting CartesianPoint to SVector

    unsorted = TypedTables.Table(detno = detno, edep = edep, pos = pos, thit = thit)
    sorting_idxs = sortperm(unsorted.detno)
    sorted = unsorted[sorting_idxs]
    grouped = TypedTables.Table(consgroupedview(sorted.detno, TypedTables.columns(sorted)))

    r_detno = similar(detno, 0)
    r_edep = similar(edep, 0)
    r_pos = similar(pos, 0)
    r_thit = similar(thit, 0)

    posunit = unit(PT)
    ustripped_cradius = ustrip(posunit, cluster_radius)

    thitunit = unit(TTT)
    ustripped_ctime = ustrip(thitunit, cluster_time)

     for d_hits_nt in grouped
        d_hits = TypedTables.Table(d_hits_nt)
        d_detno = first(d_hits.detno)
        @assert all(isequal(d_detno), d_hits.detno)

        # sort hits by time
        t_sort_idx = sortperm(d_hits.thit)
        d_hits = d_hits[t_sort_idx]
        t_vals = ustrip.(d_hits.thit)

        time_groups = Vector{Vector{Int}}()
        current_group = [1]
        start_time = t_vals[1]
        for i in 2:length(t_vals)
            if t_vals[i] - start_time ≤ ustripped_ctime
                push!(current_group, i)
            else
                current_group = finalize_group!(time_groups, current_group, i)
                start_time = t_vals[i]
            end
        end
        current_group = finalize_group!(time_groups, current_group, nothing)


        # spatial clustering within each temporal cluster
        for tg in time_groups
            t_hits = view(d_hits, tg)
            if length(t_hits) > 3
                clusters = Clustering.dbscan(hcat((ustrip.(getindex.(d_hits.pos,i)) for i in 1:3)...)', 
                ustripped_cradius, leafsize = 20, min_neighbors = 1, min_cluster_size = 1).clusters
                
                for c in clusters
                    idxs = vcat(c.boundary_indices, c.core_indices)
                    @assert length(idxs) == c.size
                    c_hits = view(d_hits, idxs)

                    push!(r_detno, d_detno)
                    esum = sum(c_hits.edep)
                    push!(r_edep, esum)
                    if esum ≈ zero(TT)
                        push!(r_pos, barycenter(c_hits.pos))
                        push!(r_thit, barycenter(c_hits.thit))
                    else
                        weights = ustrip.(Unitful.NoUnits, c_hits.edep .* inv(esum))
                        push!(r_pos, barycenter(c_hits.pos, StatsBase.Weights(weights)))
                        push!(r_thit, sum(c_hits.thit .* weights))
                    end
                end
            else
                append!(r_detno, t_hits.detno)
                append!(r_edep, t_hits.edep)
                append!(r_pos, t_hits.pos)
                append!(r_thit, t_hits.thit)
            end
        end
    end

    (detno = r_detno, edep = r_edep, pos = r_pos, thit = r_thit)
end


function cluster_detector_hits(
    table::TypedTables.Table, 
    cluster_radius::PT, 
    cluster_time::TTT
) where {PT <: RealQuantity, TTT <: RealQuantity}
    @assert :pos in TypedTables.columnnames(table) "Table has no column `pos`"
    @assert :thit in TypedTables.columnnames(table) "Table has no column `thit`"
    @assert :edep in TypedTables.columnnames(table) "Table has no column `edep`"
    @assert :detno in TypedTables.columnnames(table) "Table has no column `detno`"
    clustered_nt = map(
        evt -> cluster_detector_hits(
            evt.detno,
            evt.edep,
            [p isa CartesianPoint ? p : CartesianPoint(p...) for p in evt.pos],
            evt.thit,
            cluster_radius,
            cluster_time
        ),
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


function _group_points_by_distance(pts::AbstractVector{CartesianPoint{T}}, group_distance::T)::Tuple{Vector{Int}, Vector{Int}} where {T}
    
    n::Int = length(pts)
    
    # use BFS to find connected components
    visited        = falses(n)
    clustersidx    = similar(eachindex(pts))
    elem_ptr       = similar(eachindex(pts), n+1)
    queue          = DataStructures.CircularBuffer{Int}(n)
    
    counter        = firstindex(pts)
    Cidx           = firstindex(elem_ptr)
    elem_ptr[Cidx] = counter
    
    @inbounds for start in eachindex(pts)
        if !visited[start]
            push!(queue, start)
            visited[start] = true
            while !isempty(queue)
                node = popfirst!(queue)
                clustersidx[counter] = node
                counter += 1
                for neighbor in eachindex(pts)
                    if !visited[neighbor] && distance_squared(pts[node], pts[neighbor]) <= group_distance^2
                        push!(queue, neighbor)
                        visited[neighbor] = true
                    end
                end
            end
            Cidx += 1
            elem_ptr[Cidx] = counter
        end
    end
    clustersidx, elem_ptr[begin:Cidx]
end

function group_points_by_distance(pts::AbstractVector{CartesianPoint{T}}, group_distance::T)::Tuple{VectorOfVectors{CartesianPoint{T}}, VectorOfVectors{T}} where {T <: SSDFloat}
    clustersidx, elem_ptr = _group_points_by_distance(pts, group_distance)
    VectorOfVectors(pts[clustersidx], elem_ptr)
end

function group_points_by_distance(pts::AbstractVector{CartesianPoint{T}}, energies::AbstractVector{T}, group_distance::T)::Tuple{VectorOfVectors{CartesianPoint{T}}, VectorOfVectors{T}} where {T <: SSDFloat}
    @assert eachindex(pts) == eachindex(energies)
    clustersidx, elem_ptr = _group_points_by_distance(pts, group_distance)
    VectorOfVectors(pts[clustersidx], elem_ptr), VectorOfVectors(energies[clustersidx], elem_ptr)
end

function group_points_by_distance(pts::AbstractVector{CartesianPoint{T}}, energies::AbstractVector{T}, radius::AbstractVector{T}, group_distance::T)::Tuple{VectorOfVectors{CartesianPoint{T}}, VectorOfVectors{T}, VectorOfVectors{T}} where {T <: SSDFloat}
    @assert eachindex(pts) == eachindex(energies) == eachindex(radius)
    clustersidx, elem_ptr = _group_points_by_distance(pts, group_distance)
    VectorOfVectors(pts[clustersidx], elem_ptr), VectorOfVectors(energies[clustersidx], elem_ptr), VectorOfVectors(radius[clustersidx], elem_ptr)
end
