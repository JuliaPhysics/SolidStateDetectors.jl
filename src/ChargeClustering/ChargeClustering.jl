# This file is a part of SolidStateDetectors.jl, licensed under the MIT License (MIT).


function cluster_detector_hits(events::DetectorHitEvents, cluster_radius::RealQuantity)
    TypedTables.Table(merge(
        (evtno = events.evtno,),
        map(
            VectorOfVectors,
            Tables.columns(map(
                evt -> cluster_detector_hits(evt.detno, evt.thit, evt.edep, evt.pos, cluster_radius),
                events
            ))
        )
    ))
end


function cluster_detector_hits(
    detno::AbstractVector{<:Integer},
    thit::AbstractVector{TT},
    edep::AbstractVector{<:RealQuantity},
    pos::AbstractVector{<:StaticVector{3,PT}},
    cluster_radius::RealQuantity
) where {TT<:RealQuantity,PT <: RealQuantity}
    Table = TypedTables.Table
    unsorted = Table(detno = detno, thit = thit, edep = edep, pos = pos)
    sorting_idxs = sortperm(unsorted.detno)
    sorted = unsorted[sorting_idxs]
    grouped = Table(consgroupedview(sorted.detno, Tables.columns(sorted)))

    r_detno = similar(detno, 0)
    r_thit = similar(thit, 0)
    r_edep = similar(edep, 0)
    r_pos = similar(pos, 0)

    posunit = unit(PT)
    ustripped_cradius = ustrip(uconvert(posunit, cluster_radius))
    t_default = TT(0)

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
                
                push!(r_detno, d_detno)
                push!(r_thit, t_default)
                push!(r_edep, sum(c_hits.edep))
                push!(r_pos, mean(c_hits.pos))
            end
        else
            append!(r_detno, d_hits.detno)
            append!(r_thit, fill!(similar(d_hits.thit), t_default))
            append!(r_edep, d_hits.edep)
            append!(r_pos, d_hits.pos)
        end
    end

    (detno = r_detno, thit = r_thit, edep = r_edep, pos = r_pos)
end

export cluster_detector_hits
