# # This file is a part of SolidStateDetectors.jl, licensed under the MIT License (MIT).

function cluster_detector_hits(
        detno::AbstractVector{<:Integer},
        edep::AbstractVector{ET},
        pos::AbstractVector{<:Union{<:StaticVector{3,PT}, <:CartesianPoint{PT}}},
        thit::AbstractVector{TTT},
        cluster_radius::RealQuantity,
        cluster_time::RealQuantity = Inf*u"s"
    ) where {ET<:RealQuantity, T<:Real, PT<:Union{T,<:Unitful.Length{T}}, TT<:Real, TTT<:Union{TT,<:Unitful.Time{TT}}}

    unsorted = TypedTables.Table(detno = detno, edep = edep, pos = pos, thit = thit)
    sorting_idxs = sortperm(unsorted.detno)
    sorted = unsorted[sorting_idxs]
    grouped = TypedTables.Table(consgroupedview(sorted.detno, TypedTables.columns(sorted)))
    
    results = eltype(grouped)[]
    ustripped_cradius = _parse_value(float(T), cluster_radius, internal_length_unit)

    for d_hits_nt in grouped
        d_hits = TypedTables.Table(d_hits_nt)
        d_detno = first(d_hits.detno)
        @assert all(isequal(d_detno), d_hits.detno)

        # sort hits by time
        t_sort_idx = sortperm(d_hits.thit)
        d_hits = d_hits[t_sort_idx]
        t_vals = to_internal_units.(d_hits.thit)

        current_group = 1
        @inbounds for i in eachindex(t_vals) .+ 1
            if i > lastindex(t_vals) || t_vals[i] - t_vals[current_group] > _parse_value(float(TT), cluster_time, internal_time_unit)
                t_hits = view(d_hits, current_group:i-1)
                
                # Each time cluster gets its own result entry
                r_detno = similar(detno, 0)
                r_edep = similar(edep, 0)
                r_pos = similar(pos, 0)
                r_thit = similar(thit, 0)
                
                if length(t_hits) > 3
                    d_detno = first(t_hits.detno)
                    @assert all(isequal(d_detno), t_hits.detno)
                    clusters = Clustering.dbscan(hcat((to_internal_units.(getindex.(t_hits.pos,i)) for i in 1:3)...)', 
                        ustripped_cradius, leafsize = 20, min_neighbors = 1, min_cluster_size = 1).clusters
                    
                    for c in clusters
                        idxs = vcat(c.boundary_indices, c.core_indices)
                        @assert length(idxs) == c.size
                        c_hits = view(t_hits, idxs)
                        esum = sum(c_hits.edep)
                        push!(r_detno, d_detno)
                        push!(r_edep, esum)
                        if esum ≈ zero(ET)
                            push!(r_pos, barycenter(c_hits.pos))
                            push!(r_thit, mean(c_hits.thit))
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
                
                push!(results, (detno = r_detno, edep = r_edep, pos = r_pos, thit = r_thit))
                current_group = i
            end
        end
    end
    results
end


function cluster_detector_hits(table::TypedTables.Table, cluster_radius::RealQuantity, cluster_time::RealQuantity = Inf * u"s")
    @assert is_detector_hits_table(table) "Table does not have the correct format"
    
    # Collect all time-clustered results
    all_results = []
    evtno_col = Int[]
    
    for (evtno, evt) in enumerate(table)
        time_clusters = cluster_detector_hits(
            evt.detno,
            evt.edep,
            evt.pos,
            evt.thit,
            cluster_radius,
            cluster_time
        )
        
        append!(all_results, time_clusters)
        append!(evtno_col, fill(evtno, length(time_clusters)))
    end
    
    # Build output table
    result_columns = (
        evtno = evtno_col,
        detno = VectorOfVectors([r.detno for r in all_results]),
        thit = VectorOfVectors([r.thit for r in all_results]),
        edep = VectorOfVectors([r.edep for r in all_results]),
        pos = VectorOfVectors([r.pos for r in all_results])
    )
    
    TypedTables.Table(result_columns)
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
