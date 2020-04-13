"""
mutable struct STL{T} <: AbstractGeometry{T, 3}

Initial attempt at a complex, mesh based volume
"""

using MeshIO
using GeometryTypes
using FileIO
using NPZ
using LinearAlgebra
using Distances

struct mesh_surface{T} <: AbstractVolumePrimitive{T, 3}
    x::Array{Float32}
    y::Array{Float32}
    z::Array{Float32}
    vertices::Array{Int32}
    tri_vertices::Array{Float32}
    translate::Union{CartesianVector{T}, Missing}
end

function in(pt::CartesianPoint{T}, g::mesh_surface{T})::Bool where {T}

    use_nthreads::Int = Base.Threads.nthreads()
    
    triangles = permutedims(g.tri_vertices, [2, 3, 1])
    
    point = [pt[1] pt[2] pt[3]; pt[1] pt[2] pt[3]; pt[1] pt[2] pt[3]]

    M = triangles .- point
    M2 = (M).^2
    M3 = dropdims(sum(M2, dims=2), dims=2)
    M4 = sqrt.(M3)


    winding_sum = 0.0
    for (index_val, M4val) in enumerate(eachcol(M4)) #mapslices might work here?
        pA = det(M[:,:,index_val])
        pB = M4val[1] * M4val[2] * M4val[3]

        pC = M4val[3] * dot(M[1,:,index_val], M[2,:,index_val])
        pD = M4val[2] * dot(M[3,:,index_val], M[1,:,index_val])
        pE = M4val[1] * dot(M[2,:,index_val], M[3,:,index_val])

        pF = pC + pD + pE + pB

        winding_angle = atan(pA,pF)
        winding_sum = winding_sum + winding_angle
        
    end

    truth = ((1.9 * pi) < winding_sum < (2.1* pi))

#	open("truefalse.txt", "a") do io
#	   println(io, string(pt[1]), " ", string(pt[2]), " ", string(pt[3]), " ", string(truth))
           #write(io, pt[1], " ", pt[2], " ", pt[3], " ", truth)
#   end

    return truth

end

@inline in(pt::CylindricalPoint, g::mesh_surface)::Bool = in(CartesianPoint(pt), g)

function mesh_surface{T}(dict::Union{Dict{Any, Any}, Dict{String, Any}}, inputunit_dict::Dict{String,Unitful.Units})::mesh_surface{T} where {T <: SSDFloat}

    meshname = dict["filepath"]
    println("Reading mesh: ", meshname )
    mesh = load(meshname)
    println("STL mesh read")
    scale = [0.001, 0.001, 0.001] # should probably be user defined in a dict?
    stl_vertices = GeometryTypes.vertices(mesh)
    total_points = length(stl_vertices)
    stl_faces = GeometryTypes.faces(mesh)
    total_faces = length(stl_faces)

    stl_tri_points = zeros(Float32, total_faces, 3, 3)
    stl_face_points = zeros(Int32, total_faces, 3)
    x = zeros(Float32, total_points)
    y = zeros(Float32, total_points)
    z = zeros(Float32, total_points)
    
    index_counter = 1
    for vertex in stl_vertices
        x[index_counter] = vertex[1]*scale[1]
        y[index_counter] = vertex[2]*scale[2]
        z[index_counter] = vertex[3]*scale[3]
        index_counter = index_counter + 1
    end	

    index_counter = 1
    for stl_face in stl_faces
        stl_tri_points[index_counter, :, :] =  [stl_vertices[stl_face[1]][1]*scale[1] stl_vertices[stl_face[1]][2]*scale[2] stl_vertices[stl_face[1]][3]*scale[3]; stl_vertices[stl_face[2]][1]*scale[1] stl_vertices[stl_face[2]][2]*scale[2] stl_vertices[stl_face[2]][3]*scale[3]; stl_vertices[stl_face[3]][1]*scale[1] stl_vertices[stl_face[3]][2]*scale[2] stl_vertices[stl_face[3]][3]*scale[3]]
        stl_face_points[index_counter, :] =  [stl_face[1], stl_face[2], stl_face[3]]
        
        index_counter = index_counter + 1

    end 

    if haskey(dict, "translate")
        translate = CartesianVector{T}(
            haskey(dict["translate"],"x") ? geom_round(ustrip(uconvert(u"m", T(dict["translate"]["x"]) * inputunit_dict["length"] ))) : 0.0,
            haskey(dict["translate"],"y") ? geom_round(ustrip(uconvert(u"m", T(dict["translate"]["y"]) * inputunit_dict["length"] ))) : 0.0,
            haskey(dict["translate"],"z") ? geom_round(ustrip(uconvert(u"m", T(dict["translate"]["z"]) * inputunit_dict["length"] ))) : 0.0)
    else
        translate = CartesianVector{T}(0.0, 0.0, 0.0)
    end

    return mesh_surface{T}(
        x,
        y,
        z,
        stl_face_points,
        stl_tri_points,
        translate)
end

function Geometry(T::DataType, t::Val{:mesh_surface}, dict::Dict{Any, Any}, inputunit_dict::Dict{String,Unitful.Units})
    return mesh_surface{T}(dict, inputunit_dict)
end

function get_important_points(g::mesh_surface{T})::NTuple{3, Vector{T}} where {T <: SSDFloat}
    v1::Vector{T} = T[g.x[1], g.x[2]] #[g.x[1], g.x[2]]
    v2::Vector{T} = T[g.y[1], g.y[2]]
    v3::Vector{T} = T[g.z[1], g.z[2]]
    return v1, v2, v3
end

function get_important_points(g::mesh_surface{T}, ::Val{:r})::Vector{T} where {T <: SSDFloat}
    return geom_round.(T.(abs.([g.x[1], g.x[2], g.y[1], g.y[2]])))
end
function get_important_points(g::mesh_surface{T}, ::Val{:φ})::Vector{T} where {T <: SSDFloat}
    return T[]
end
function get_important_points(g::mesh_surface{T}, ::Val{:z})::Vector{T} where {T <: SSDFloat}
    return geom_round.(T[g.z[1], g.z[2]])
end
function get_important_points(g::mesh_surface{T}, ::Val{:x})::Vector{T} where {T <: SSDFloat}
    return geom_round.(T[g.x[1], g.x[2]])
end
function get_important_points(g::mesh_surface{T}, ::Val{:y})::Vector{T} where {T <: SSDFloat}
    return geom_round.(T[g.y[1], g.y[2]])
end


function vertices(cb::mesh_surface{T})::Vector{CartesianPoint{T}} where {T <: SSDFloat}

    v = CartesianPoint{T}[]

    #println("Making vertices")
    index_counter = 1
    for x_val in cb.x
        push!(v, CartesianPoint{T}(cb.x[index_counter], cb.y[index_counter], cb.z[index_counter]))

        index_counter = index_counter + 1
    end	

    ismissing(cb.translate) ? nothing : v = map(x -> x + cb.translate, v)
    #println("Vertices produced")
    return v
end

function LineSegments(cb::mesh_surface{T})::Vector{LineSegment{T, 3, :cartesian}} where {T}
    v::Vector{CartesianPoint{T}} = vertices(cb)

    ls = LineSegment{T, 3, :cartesian}[]

    for face_indices in eachrow(cb.vertices)
        push!(ls, LineSegment(v[face_indices[1]], v[face_indices[2]]))
        push!(ls, LineSegment(v[face_indices[2]], v[face_indices[3]]))
        push!(ls, LineSegment(v[face_indices[3]], v[face_indices[1]]))
    end	

    return ls
end

@recipe function f(cb::mesh_surface{T}) where {T <: SSDFloat}
    label-->"mesh_surface"
    ls = LineSegments(cb)
    @series begin
        ls
    end
end

function sample(cb::mesh_surface{T}, stepsize::Vector{T})  where {T <: SSDFloat}
    samples  = CartesianPoint{T}[]
    for x in cb.x[1] : stepsize[1] : cb.x[2]
        for y in cb.y[1] :stepsize[2] : cb.y[2]
            for z in cb.z[1] : stepsize[3] : cb.z[2]
                push!(samples,CartesianPoint{T}(x, y, z))
            end
        end
    end
    ismissing(cb.translate) ? nothing : samples = map(x -> x + cb.translate, samples)
    return samples
end

function (+)(b::mesh_surface{T}, translate::Union{CartesianVector{T},Missing})::mesh_surface{T} where {T <: SSDFloat}
    if ismissing(translate)
        return b
    elseif ismissing(b.translate)
        return mesh_surface(b.x, b.y, b.z, b.vertices, translate)
    else
        return mesh_surface(b.x, b.y, b.z, b.vertices, b.tri_vertices, b.translate + translate)
    end
 end