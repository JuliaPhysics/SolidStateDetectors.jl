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

"""
    nidx( rbidx::Int, ::Val{true}, ::Val{true})::Int

first type argument:  type of the original point (if even points -> `Val{true}()`, else `Val{false}()`)
second type argument: is sum of other point indices even or odd -> (if sum is even -> `Val{true}()`, else `Val{false}()`)
"""
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

"""
    get_rbidx_right_neighbour(rbidx::Int, ::Val{true}, ::Val{true})::Int

needs docu...
"""
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


"""
    RBExtBy2Array( et::Type, g::Grid{T, N, :Cylindrical} )::Array{et, N + 1} where {T, N}

Returns a RedBlack array for the grid `g`. 
"""
function RBArray( et::Type, g::Grid{T, N, :Cylindrical} )::Array{et, N + 1} where {T, N}
    nr, nθ, nz = size(g)
    # new ordering in memory: r, θ, z -> z, θ, r (so inner loop goes over z)
    return zeros(T, div(nz, 2) + mod(nz, 2), nθ, nr, 2)
end
function RBArray( a::Array{T, N}, grid::Grid{TG, N, :Cylindrical} )::Array{T, N + 1} where {T, N, TG}
    rbarray::Array{T, N + 1} = RBArray(T, grid)
    for iz in axes(a, 3)
        irbz::Int = div(iz, 2) + mod(iz, 2) 
        for iθ in axes(a, 2)
            idxsum::Int = iz + iθ
            for ir in axes(a, 1)
                rbi::Int = iseven(idxsum + ir) ? rb_even::Int : rb_odd::Int
                rbarray[ irbz, iθ, ir, rbi ] = a[ir, iθ, iz]
            end
        end
    end
    return rbarray
end

# """
#     RBExtBy1Array( et::Type, g::Grid{T, N, :Cylindrical} )::Array{et, N + 1} where {T, N}

# Returns a RedBlack array for the grid `g`. The RedBlack array is extended in its size by 1 in each geometrical dimension.
# """
# function RBExtBy1Array( et::Type, g::Grid{T, N, :Cylindrical} )::Array{et, N + 1} where {T, N}
#     nr, nθ, nz = size(g)
#     # new ordering in memory: r, θ, z -> z, θ, r (so inner loop goes over z)
#     return zeros(T, div(nz, 2) + mod(nz, 2) + 1, nθ + 1, nr + 1, 2)
# end
# function RBExtBy1Array( a::Array{T, N}, grid::Grid{TG, N, :Cylindrical} )::Array{T, N + 1} where {T, N, TG}
#     rbarray::Array{T, N + 1} = RBExtBy1Array(T, grid)
#     for iz in axes(a, 3)
#         irbz::Int = div(iz, 2) + mod(iz, 2)  
#         for iθ in axes(a, 2)
#             irbθ::Int = iθ
#             idxsum::Int = iz + iθ
#             for ir in axes(a, 1)
#                 irbr::Int = ir 
#                 rbi::Int = iseven(idxsum + ir) ? rb_even::Int : rb_odd::Int
#                 rbarray[ irbz, irbθ, irbr, rbi ] = a[ir, iθ, iz]
#             end
#         end
#     end
#     return rbarray
# end

"""
    RBExtBy2Array( et::Type, g::Grid{T, N, :Cylindrical} )::Array{et, N + 1} where {T, N}

Returns a RedBlack array for the grid `g`. The RedBlack array is extended in its size by 2 in each geometrical dimension.
"""
function RBExtBy2Array( et::Type, g::Grid{T, N, :Cylindrical} )::Array{et, N + 1} where {T, N}
    nr, nθ, nz = size(g)
    # new ordering in memory: r, θ, z -> z, θ, r (so inner loop goes over z)
    return zeros(T, div(nz, 2) + mod(nz, 2) + 2, nθ + 2, nr + 2, 2)
end
function RBExtBy2Array( a::Array{T, N}, grid::Grid{TG, N, :Cylindrical} )::Array{T, N + 1} where {T, N, TG}
    rbarray::Array{T, N + 1} = RBExtBy2Array(T, grid)
    for iz in axes(a, 3)
        irbz::Int = rbidx(iz) 
        for iθ in axes(a, 2)
            irbθ::Int = iθ + 1
            idxsum::Int = iz + iθ
            for ir in axes(a, 1)
                irbr::Int = ir + 1
                rbi::Int = iseven(idxsum + ir) ? rb_even::Int : rb_odd::Int
                rbarray[ irbz, irbθ, irbr, rbi ] = a[ir, iθ, iz]
            end
        end
    end
    return rbarray
end

# function NArray( rbarray::Array{T, N1}, grid::Grid{TG, N2, :Cylindrical} )::Array{T, N2} where {T, TG, N1, N2}
#     a::Array{T, N2} = Array{T, N2}(undef, size(grid)...)
#     if size(rbarray, 3) == size(grid, 1) # RBArray (not extended)
#         for iz in axes(a, 3)
#             irbz::Int = div(iz, 2) + mod(iz, 2) 
#             for iθ in axes(a, 2)
#                 idxsum::Int = iz + iθ
#                 for ir in axes(a, 1)
#                     rbi::Int = iseven(idxsum + ir) ? rb_even::Int : rb_odd::Int
#                     a[ir, iθ, iz] = rbarray[ irbz, iθ, ir, rbi ]
#                 end
#             end
#         end
#     elseif size(rbarray, 3) == size(grid, 1) + 2 # extended RBArray by 2
#         for iz in axes(a, 3)
#             irbz::Int = rbidx(iz) 
#             for iθ in axes(a, 2)
#                 irbθ::Int = iθ + 1
#                 idxsum::Int = iz + iθ
#                 for ir in axes(a, 1)
#                     irbr::Int = ir + 1
#                     rbi::Int = iseven(idxsum + ir) ? rb_even::Int : rb_odd::Int
#                     a[ir, iθ, iz] = rbarray[ irbz, irbθ, irbr, rbi ]
#                 end
#             end
#         end
#     else
#         error("`NArray()` not defined for this array size.")
#     end
#     return a
# end

