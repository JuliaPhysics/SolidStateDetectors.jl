# # This file is a part of SolidStateDetectors.jl, licensed under the MIT License (MIT).

function cluster_detector_hits(
        detno::AbstractVector{<:Integer},
        edep::AbstractVector{ET},
        pos::AbstractVector{<:Union{<:StaticVector{3,PT}, <:CartesianPoint{PT}}},
        thit::AbstractVector{TTT},
        cluster_radius::RealQuantity
    ) where {ET<:RealQuantity, T<:Real, PT<:Union{T,<:Unitful.Length{T}}, TT<:Real, TTT<:Union{TT,<:Unitful.Time{TT}}}

    unsorted = TypedTables.Table(detno = detno, edep = edep, pos = pos, thit = thit)
    sorting_idxs = sortperm(unsorted.detno)
    sorted = unsorted[sorting_idxs]
    grouped = TypedTables.Table(consgroupedview(sorted.detno, TypedTables.columns(sorted)))

    r_detno = similar(detno, 0)
    r_edep = similar(edep, 0)
    r_pos = similar(pos, 0)
    r_thit = similar(thit, 0)

    ustripped_cradius = _parse_value(float(T), cluster_radius, internal_length_unit)
    
    for d_hits_nt in grouped
        d_hits = TypedTables.Table(d_hits_nt)
        d_detno = first(d_hits.detno)
        @assert all(isequal(d_detno), d_hits.detno)
        if length(d_hits) > 3
            
            clusters = Clustering.dbscan(hcat((to_internal_units.(getindex.(d_hits.pos,i)) for i in 1:3)...)', 
                ustripped_cradius, leafsize = 20, min_neighbors = 1, min_cluster_size = 1).clusters
            for c in clusters
                idxs = vcat(c.boundary_indices, c.core_indices)
                @assert length(idxs) == c.size
                c_hits = view(d_hits, idxs)
                
                push!(r_detno, d_detno)
                esum = sum(c_hits.edep)
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
            append!(r_detno, d_hits.detno)
            append!(r_edep, d_hits.edep)
            append!(r_pos, d_hits.pos)
            append!(r_thit, d_hits.thit)
        end
    end

    (detno = r_detno, edep = r_edep, pos = r_pos, thit = r_thit)
end


function cluster_detector_hits(table::TypedTables.Table, cluster_radius::RealQuantity, cluster_time::RealQuantity = Inf * u"s")
    @assert is_detector_hits_table(table) "Table does not have the correct format"

    if isfinite(cluster_time)
        table = split_by_time(table, cluster_time)
    end
    
    clustered_nt = map(
        evt -> cluster_detector_hits(evt.detno, evt.edep, evt.pos, evt.thit, cluster_radius),
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

function group_points_by_distance(pts::AbstractVector{CartesianPoint{T}}, group_distance::T)::VectorOfVectors{CartesianPoint{T}} where {T <: SSDFloat}
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
