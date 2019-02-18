
struct Cylinder{T} <: AbstractGeometry{T, 3}
    org::CartesianPoint{T}
    ext::CartesianPoint{T} # optional, but saved computation time
    dir::CartesianVector{T} # not normalized -> norm(dir) = length
    radius::T
    length::T # optional, but saved computation time 
end

function Cylinder{T}(org::CartesianPoint{T}, dir::CartesianVector{T}, r::T) where {T}
    Cylinder{T}(org, org + dir, dir, r, norm(dir))
end

function in(pt::CartesianPoint{T}, cyl::Cylinder{T})::Bool where {T <: Real}
    shift::CartesianVector{T} = pt - cyl.org
    return  (shift ⋅ cyl.dir >= 0) && 
            ((pt - cyl.ext) ⋅ cyl.dir <= 0) &&
            (norm(cross(shift, cyl.dir)) / cyl.length <= cyl.radius )
end
