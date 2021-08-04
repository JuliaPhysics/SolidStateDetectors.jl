"""
    struct NBodyChargeCloud{T} <: AbstractChargeCloud

Struct which defines a charge cloud consisting of multiple point-like charge carriers,
initially distributed around a given center.

## Fields
* `points::Vector{CartesianPoint{T}}`: Positions of the charge carriers that are part of the charge cloud.
* `energies::Vector{T}`: Energies of the respective charge carriers, in the same order as `points`.
* `shell_structure::Vector{Type{<:AbstractChargeCloud}}`: Initial geometry of the charge carriers around the
    `center` point, relevant for plotting.

See also [`create_charge_cloud`](@ref).
"""
struct NBodyChargeCloud{T} <: AbstractChargeCloud
    points::Vector{CartesianPoint{T}}
    energies::Vector{T}
    shell_structure::Vector{Type{<:AbstractChargeCloud}}
end

@recipe function f(nbc::NBodyChargeCloud{T}; connect = true, markersize = 10) where {T}
    
    seriescolor --> :blue
    
    @series begin
        seriescolor --> seriescolor
        connect --> connect
        seriestype --> :scatter
        markersize --> markersize
        markerstrokewidth --> 0
        label --> "NBodyChargeCloud"
        []
    end
    
    points = nbc.points
    vertex_no = 0
    for (shell, shell_structure) in enumerate(nbc.shell_structure)
        vertices = get_vertices(shell_structure)
        @series begin
            markersize --> markersize * exp(-(shell - 1))
            label --> ""
            seriescolor --> seriescolor
            shell_structure(points[vertex_no+1:vertex_no+vertices])
        end
        vertex_no += vertices
    end
end
