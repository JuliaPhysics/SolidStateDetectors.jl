# # This file is a part of SolidStateDetectors.jl, licensed under the MIT License (MIT).

# breaking changes in Clustering v0.14 -> v0.15
@inline _get_clusters(clusters::Clustering.DbscanResult) = clusters.clusters
@inline _get_clusters(clusters::Vector{Clustering.DbscanCluster}) = clusters

function cluster_detector_hits(
    detno::AbstractVector{<:Integer},
    edep::AbstractVector{TT},
    pos::AbstractVector{<:StaticVector{3,PT}},
    cluster_radius::RealQuantity
) where {TT<:RealQuantity, PT <: RealQuantity}
    Table = TypedTables.Table
    unsorted = Table(detno = detno, edep = edep, pos = pos)
    sorting_idxs = sortperm(unsorted.detno)
    sorted = unsorted[sorting_idxs]
    grouped = Table(consgroupedview(sorted.detno, TypedTables.columns(sorted)))

    r_detno = similar(detno, 0)
    r_edep = similar(edep, 0)
    r_pos = similar(pos, 0)

    posunit = unit(PT)
    ustripped_cradius = ustrip(posunit, cluster_radius)
    
    for d_hits_nt in grouped
        d_hits = Table(d_hits_nt)
        d_detno = first(d_hits.detno)
        @assert all(isequal(d_detno), d_hits.detno)
        if length(d_hits) > 3
            clusters = Clustering.dbscan(ustrip.(flatview(d_hits.pos)), ustripped_cradius, leafsize = 20, min_neighbors = 1, min_cluster_size = 1)
            for c in _get_clusters(clusters)
                idxs = vcat(c.boundary_indices, c.core_indices)
                @assert length(idxs) == c.size
                c_hits = view(d_hits, idxs)
                
                push!(r_detno, d_detno)
                esum_u = sum(c_hits.edep)
                push!(r_edep, esum_u)
                esum = ustrip(esum_u)
                if esum ≈ 0
                    push!(r_pos, mean(c_hits.pos))
                else
                    weights = ustrip.(c_hits.edep) .* inv(esum)
                    push!(r_pos, sum(c_hits.pos .* weights))
                end
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
