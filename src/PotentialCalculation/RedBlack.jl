abstract type AbstractRBArray{T, N} <: AbstractArray{T, N} end

abstract type RBEven end
abstract type RBOdd end

const rb_even = Int(2)
const rb_odd  = Int(1)

const rb_bool_even = true
const rb_bool_odd  = false


@inline function rbidx( nidx::Int )::Int
    return div(nidx, 2) + mod(nidx, 2) + 1
end

# """
#     nidx( rbidx::Int, ::Val{true}, ::Val{true})::Int
# 
# first type argument:  type of the orgal point (for even points -> `Val{true}()`, else `Val{false}()`)
# second type argument: is sum of other point indices even or odd -> (if sum is even -> `Val{true}()`, else `Val{false}()`)
# """
@inline function nidx( rbidx::Int, ::Val{true}, ::Val{true})::Int
   return (rbidx - 1) * 2 
end
@inline function nidx( rbidx::Int, ::Val{true}, ::Val{false})::Int
   return (rbidx - 1) * 2 - 1
end
@inline function nidx( rbidx::Int, ::Val{false}, ::Val{true})::Int
   return (rbidx - 1) * 2 - 1
end
@inline function nidx( rbidx::Int, ::Val{false}, ::Val{false})::Int
   return (rbidx - 1) * 2 
end
@inline function nidx( rbidx::Int, b1::Bool, b2::Bool)::Int
    if xor(b1, b2)
        (rbidx - 1) * 2 - 1
    else
        (rbidx - 1) * 2 
    end
end

# """
#     get_rbidx_right_neighbour(rbidx::Int, ::Val{true}, ::Val{true})::Int
# 
# needs docu...
# """
@inline function get_rbidx_right_neighbour(rbidx::Int, ::Val{true}, ::Val{true})::Int
    return rbidx + 1
end
@inline function get_rbidx_right_neighbour(rbidx::Int, ::Val{true}, ::Val{false})::Int
    return rbidx 
end
@inline function get_rbidx_right_neighbour(rbidx::Int, ::Val{false}, ::Val{true})::Int
    return rbidx 
end
@inline function get_rbidx_right_neighbour(rbidx::Int, ::Val{false}, ::Val{false})::Int
    return rbidx + 1
end
@inline function get_rbidx_right_neighbour(rbidx::Int, b1::Bool, b2::Bool)::Int
    if xor(b1, b2)
        rbidx
    else
        rbidx + 1
    end
end


# """
#     RBExtBy2Array( T::Type, grid::Grid{TG, N, Cylindrical} )::Array{T, N + 1} where {TG, N}
# 
# Returns a RedBlack array for the `grid`. 
# """
function RBArray( T::Type, grid::Grid{TG, N, Cylindrical} )::Array{T, N + 1} where {TG, N}
    nr, nφ, nz = size(grid)
    # new ordering in memory: r, φ, z -> z, φ, r (so inner loop goes over z)
    return zeros(T, div(nz, 2) + mod(nz, 2), nφ, nr, 2)
end
function RBArray( a::Array{T, N}, grid::Grid{TG, N, Cylindrical} )::Array{T, N + 1} where {T, N, TG}
    rbarray::Array{T, N + 1} = RBArray(T, grid)
    for iz in axes(a, 3)
        irbz::Int = div(iz, 2) + mod(iz, 2) 
        for iφ in axes(a, 2)
            idxsum::Int = iz + iφ
            for ir in axes(a, 1)
                rbi::Int = iseven(idxsum + ir) ? rb_even::Int : rb_odd::Int
                rbarray[ irbz, iφ, ir, rbi ] = a[ir, iφ, iz]
            end
        end
    end
    return rbarray
end

# """
#     RBExtBy2Array( T::Type, grid::Grid{TG, N, Cylindrical} )::Array{T, N + 1} where {TG, N}
# 
# Returns a RedBlack array for the `grid`. The RedBlack array is extended in its size by 2 in each geometrical dimension.
# """
function RBExtBy2Array( T::Type, grid::Grid{TG, N, Cylindrical} )::Array{T, N + 1} where {TG, N}
    nr, nφ, nz = size(grid)
    # new ordering in memory: r, φ, z -> z, φ, r (so inner loop goes over z)
    return zeros(T, div(nz, 2) + mod(nz, 2) + 2, nφ + 2, nr + 2, 2)
end
function RBExtBy2Array( a::Array{T, N}, grid::Grid{TG, N, Cylindrical} )::Array{T, N + 1} where {T, N, TG}
    rbarray::Array{T, N + 1} = RBExtBy2Array(T, grid)
    for iz in axes(a, 3)
        irbz::Int = rbidx(iz) 
        for iφ in axes(a, 2)
            irbφ::Int = iφ + 1
            idxsum::Int = iz + iφ
            for ir in axes(a, 1)
                irbr::Int = ir + 1
                rbi::Int = iseven(idxsum + ir) ? rb_even::Int : rb_odd::Int
                rbarray[ irbz, irbφ, irbr, rbi ] = a[ir, iφ, iz]
            end
        end
    end
    return rbarray
end


# """
#     RBExtBy2Array( T::Type, grid::CartesianGrid3D{TG} )::Array{T, 4} where {TG}
# 
# Returns a RedBlack array for the `grid`. The RedBlack array is extended in its size by 2 in each geometrical dimension.
# """
function RBExtBy2Array( T::Type, grid::CartesianGrid3D{TG} )::Array{T, 4} where {TG}
    nx, ny, nz = size(grid)
    return zeros(T, div(nx, 2) + mod(nx, 2) + 2, ny + 2, nz + 2, 2)
end
function RBExtBy2Array( a::Array{T, 3}, grid::Grid{TG, 3, Cartesian} )::Array{T, 4} where {T, TG}
    rbarray::Array{T, 4} = RBExtBy2Array(T, grid)
    for iz in axes(a, 3)
        irbz::Int = iz + 1
        for iy in axes(a, 2)
            irby::Int = iy + 1
            idxsum::Int = iz + iy
            for ix in axes(a, 1)
                irbx::Int = rbidx(ix) 
                rbi::Int = iseven(idxsum + ix) ? rb_even::Int : rb_odd::Int
                rbarray[ irbx, irby, irbz, rbi ] = a[ix, iy, iz]
            end
        end
    end
    return rbarray
end


