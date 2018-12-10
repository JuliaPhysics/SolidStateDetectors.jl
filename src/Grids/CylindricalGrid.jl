"""
    CylindricalGrid{T<:AbstractFloat} <: Grid

Stores the electric potential on a three dimensional cylindrical grid and also the position of the grid points.

# Fields
- `cyclic::T`: Stores the periodicity of the system. E.g. 2π.
- `r::T`: The r-axis of the grid in `m`.
- `φ::T`: The φ-axis of the grid in `rad`.
- `z::T`: The r-axis of the grid in `m`.
- `potential::Array{T, 3}`: The potential values in `V` of each grid point.
"""
mutable struct CylindricalGrid{T<:AbstractFloat} <: Grid
    cyclic::T
    "Units: m"
    r::Array{T, 1}
    "Units: rad"
    φ::Array{T, 1}  # use radiance
    "Units: m"
    z::Array{T, 1}
    "Units: V"
    potential::Array{T, 3}

    function CylindricalGrid{T}(cyclic::Real, r::AbstractArray, φ::AbstractArray, z::AbstractArray; only_2d::Bool = false) where {T<:AbstractFloat}
        if !only_2d @assert iseven(length(φ)) "Unhandled boundary conditions, due to uneven grid in φ." end
        return  new{T}( cyclic, r, φ, z, zeros(T, length(r), length(φ), length(z)) )
    end
end


function CylindricalGrid(nt::NamedTuple)
    T = typeof(ustrip(nt.values[1]))
    grid = CylindricalGrid{T}(
        T(ustrip(uconvert(u"rad", nt.periodicity.phi))),
        convert(Array{T}, ustrip.(uconvert.(u"m", nt.edges.r))),
        convert(Array{T}, ustrip.(uconvert.(u"rad", nt.edges.phi))),
        convert(Array{T}, ustrip.(uconvert.(u"m", nt.edges.z))),
        only_2d = length(nt.edges.phi) == 1,
    )
    grid.potential = convert(Array{T}, ustrip.(uconvert.(u"V", nt.values)))
    grid
end

Base.convert(T::Type{CylindricalGrid}, x::NamedTuple) = T(x)


function NamedTuple(grid::CylindricalGrid)
    T = eltype(grid.potential)
    (
        values = grid.potential * u"V",
        edges = (
            r = grid.r * u"m",
            phi = grid.φ * u"rad",
            z = grid.z * u"m"
        ),
        periodicity = (
            phi = grid.cyclic * u"rad",
        )
    )
end

Base.convert(T::Type{NamedTuple}, x::CylindricalGrid) = T(x)


const ElectricPotential = CylindricalGrid
const WeightingPotential = CylindricalGrid

function Base.println(io::IO, grid::CylindricalGrid)
    println("CylindricalGrid:")
    println("cyclic: $(grid.cyclic)°")
    println("potential: $(typeof(grid.potential)) - dims: $(size(grid.potential))")
    if length(grid.r) > 0
        println("r: $(minimum(grid.r)) - $(maximum(grid.r)), $(length(grid.r)) points")
    else
        println("r: empty")
    end
    if length(grid.φ) > 0
        println("φ: $(minimum(grid.φ)) - $(maximum(grid.φ)), $(length(grid.φ)) points")
    else
        println("φ: empty")
    end
    if length(grid.z) > 0
        println("z: $(minimum(grid.z)) - $(maximum(grid.z)), $(length(grid.z)) points")
    else
        println("z: empty")
    end
end
function Base.print(io::IO, grid::CylindricalGrid) println(io, grid) end


function Base.show(io::IO, grid::CylindricalGrid) println(io, grid) end

function Base.show(io::IO, ::MIME"text/plain", grid::CylindricalGrid)
    show(io, grid)
end
function Base.display(io::IO, grid::CylindricalGrid) println(io, grid) end





function CylindricalGrid(detector::SolidStateDetector, r_arr::AbstractArray, φ_arr::AbstractArray, z_arr::AbstractArray, cyclic::Real; only_2d::Bool = false)::CylindricalGrid
    T = get_precision_type(detector)

    important_r_points = sort(round.(get_important_r_points(detector), sigdigits=6))
    important_φ_points = T[]#!only_2d ? sort(get_important_φ_points(detector)) : T[]
    important_z_points = sort(round.(get_important_z_points(detector), sigdigits=6)) #T[]

    for r in important_r_points
        if !in(r, r_arr) push!(r_arr, r) end
    end
    for φ in important_φ_points
        if 0 <= φ < cyclic
            already_inside = false
            for aφ in φ_arr
                if (round(aφ, sigdigits=6) ≈ round(φ, sigdigits=6))
                    already_inside = true
                end
            end
            if !already_inside
                push!(φ_arr, φ)
            end
        end
    end
    for z in important_z_points
        if !in(z, z_arr) push!(z_arr, z) end
    end

    sort!(r_arr)
    sort!(φ_arr)
    sort!(z_arr)

    if isodd(length(φ_arr)) && (!only_2d)
        max_diff_idx = findmax(diff(φ_arr))[2]
        push!(φ_arr, 0.5 * (φ_arr[max_diff_idx + 1] + φ_arr[max_diff_idx]))
        sort!(φ_arr)
    end
    if isodd(length(z_arr))
        max_diff_idx = findmax(diff(z_arr))[2]
        push!(z_arr, 0.5 * (z_arr[max_diff_idx + 1] + z_arr[max_diff_idx]))
        sort!(z_arr)
    end

    return CylindricalGrid{T}(cyclic, r_arr, φ_arr, z_arr, only_2d=only_2d)
end

function CylindricalGrid(detector::SolidStateDetector; init_grid_spacing::Real=2e-3, cyclic::Real=2π, only_2d::Bool = false)::CylindricalGrid
    T = get_precision_type(detector)

    only_2d::Bool = cyclic == 0 ? true : false

    grid_spacing::T = init_grid_spacing

    r_start::T = 0
    r_stop::T = detector.crystal_radius + 4 * grid_spacing
    r_step::T = grid_spacing # 1mm
    r_length::Int = div(r_stop - r_start, r_step)
    r_arr::Vector{T} = collect(range(r_start, stop=r_stop, length=r_length))

    φ_start::T = 0
    φ_step::T = grid_spacing / r_arr[2]
    
    # φ_length::Int = div(φ_stop - φ_start, φ_step)
    φ_length = 8
    φ_step = cyclic / φ_length 
    φ_stop::T = cyclic - φ_step

    if isodd(φ_length) φ_length += 1 end
    if φ_length == 0 φ_length = 2 end
    φ_arr::Vector{T} = if only_2d
        T[0]
    else
        collect(range(φ_start, stop=φ_stop, length=φ_length) )
    end

    z_start::T =  -4 * grid_spacing
    z_step::T = grid_spacing # 1 mm
    z_stop::T = detector.crystal_length + 4 * grid_spacing
    z_length::Int = div(z_stop - z_start, z_step)
    if isodd(z_length)
        z_length + 1
        z_stop += z_step
    end

    z_arr::Vector{T} = collect(range( z_start, stop=z_stop, length=z_length) )

    return CylindricalGrid(detector, r_arr, φ_arr, z_arr, cyclic, only_2d=only_2d)
end


function add_points_in_z_and_interpolate!(grid::CylindricalGrid, idcs::AbstractArray)::Nothing
    T = eltype(grid.potential)
    n = length(idcs)
    idcs = collect(idcs)
    if n == 0 return nothing end
    if !issorted(idcs) idcs = sort(idcs) end
    if idcs[end] >= length(grid.z)
        error("Can not append point in z at the end. Only in between.")
    end

    ranges = [ 1:idcs[1] ]
    for i in (2:n)
        push!(ranges, idcs[i-1]+1:idcs[i])
    end
    push!(ranges, idcs[end]+1:length(grid.z))

    grid.z = cat(grid.z, zeros(T, n), dims=1)

    for (i, range) in enumerate(reverse(ranges))
        grid.z[range .+ (n+1-i)] = grid.z[range]
    end
    for (i, range) in enumerate(reverse(ranges[2:end]))
        grid.z[first(range .+ (n+1-i))-1] = 0.5*(grid.z[first(range .+ (n+1-i))]+grid.z[first(range .+ (n+1-i))-2])
    end

    # 1: extend arrays
     grid.potential = cat( grid.potential, zeros(T, size( grid.potential, 1), size( grid.potential, 2), n), dims=3)
    # 2: copy old values to new locations
    for (i, range) in enumerate(reverse(ranges))
        grid.potential[:, :, range .+ (n+1-i)] = grid.potential[:, :, range]
    end
     # 3: interpolate new positions
    for (i, range) in enumerate(reverse(ranges[2:end]))
        grid.potential[:, :, first(range .+ (n+1-i))-1] = 0.5 * ( grid.potential[:, :, first(range .+ (n+1-i))] + grid.potential[:, :, first(range .+ (n+1-i))-2])
    end

    return nothing
end

function add_points_in_r_and_interpolate!(grid::CylindricalGrid, idcs::AbstractArray)::Nothing
    T = eltype(grid.potential)
    n = length(idcs)
    idcs = collect(idcs)
    if n == 0 return nothing end
    if !issorted(idcs) idcs = sort(idcs) end
    if idcs[end] >= length(grid.r)
        error("Can not append point in r at the end. Only in between.")
    end

    ranges = [ 1:idcs[1] ]
    for i in (2:n)
        push!(ranges, idcs[i-1]+1:idcs[i])
    end
    push!(ranges, idcs[end]+1:length(grid.r))

    grid.r = cat(grid.r, zeros(T, n), dims=1)

    for (i, range) in enumerate(reverse(ranges))
        grid.r[range .+ (n+1-i)] = grid.r[range]
    end
    for (i, range) in enumerate(reverse(ranges[2:end]))
        grid.r[first(range .+ (n+1-i))-1] = 0.5*(grid.r[first(range .+ (n+1-i))]+grid.r[first(range .+ (n+1-i))-2])
    end

    # 1: extend arrays
    grid.potential = cat( grid.potential, zeros(T, n, size( grid.potential, 2), size( grid.potential, 3)), dims=1)
    # 2: copy old values to new locations
    for (i, range) in enumerate(reverse(ranges))
        grid.potential[range .+ (n+1-i), :, :] = grid.potential[range, :, :]
    end
     # 3: interpolate new positions
    for (i, range) in enumerate(reverse(ranges[2:end]))
        grid.potential[first(range .+ (n+1-i))-1, :, :] = 0.5 * ( grid.potential[first(range .+ (n+1-i)), :, :] + grid.potential[first(range .+ (n+1-i))-2, :, :])
    end

    return nothing
end

function add_points_in_φ_and_interpolate!(grid::CylindricalGrid, idcs::AbstractArray)::Nothing
    T = eltype(grid.potential)
    n = length(idcs)
    idcs = collect(idcs)
    if n == 0 return nothing end
    if !issorted(idcs) idcs = sort(idcs) end

    if idcs[end] > length(grid.φ)
        error("Can not append array in φ at index $(idcs[end]). (Cylcic: $(grid.cyclic)°)")
    end

    append_after_last_index = false
    if idcs[end] == length(grid.φ)
        append_after_last_index = true
        idcs = idcs[1:end-1]
        n -= 1
    end

    ranges = []
    if length(idcs) > 0
        ranges = [ 1:idcs[1] ]
        for i in (2:n)
            push!(ranges, idcs[i-1]+1:idcs[i])
        end
        if idcs[end] < length(grid.φ)
            push!(ranges, idcs[end]+1:length(grid.φ))
        end
    end

    if append_after_last_index
        grid.φ = cat( grid.φ, zeros(T, n+1), dims=1)
    else
        grid.φ = cat( grid.φ, zeros(T, n), dims=1)
    end

    for (i, range) in enumerate(reverse(ranges))
        grid.φ[range .+ (n+1-i)] = grid.φ[range]
    end
    for (i, range) in enumerate(reverse(ranges[2:end]))
        grid.φ[first(range .+ (n+1-i))-1] = 0.5*(grid.φ[first(range .+ (n+1-i))]+grid.φ[first(range .+ (n+1-i))-2])
    end

    if append_after_last_index
        grid.φ[end] = 0.5 * (grid.cyclic + grid.φ[end-1])
    end

    # 1: extend arrays
    if append_after_last_index
        grid.potential = cat( grid.potential, zeros(T, size( grid.potential, 1), n + 1, size( grid.potential, 3)), dims=2 )
    else
        grid.potential = cat( grid.potential, zeros(T, size( grid.potential, 1), n, size( grid.potential, 3)), dims=2)
    end
    # 2: copy old values to new locations
    for (i, range) in enumerate(reverse(ranges))
        grid.potential[:, range .+ (n+1-i), :] = grid.potential[:, range, :]
    end
     # 3: interpolate new positions
    for (i, range) in enumerate(reverse(ranges[2:end]))
        grid.potential[: ,first(range .+ (n+1-i))-1, :] = 0.5 * ( grid.potential[:, first(range .+ (n+1-i)), :] + grid.potential[:, first(range .+ (n+1-i))-2, :])
    end
    if append_after_last_index
        grid.potential[: ,end, :] = 0.5 * ( grid.potential[:, end-1, :] + grid.potential[:, 1, :])
    end

    return nothing
end


function add_points_and_interpolate!(grid::CylindricalGrid, idcs_r::AbstractArray, idcs_φ::AbstractArray, idcs_z::AbstractArray)::Nothing
    add_points_in_r_and_interpolate!(grid, idcs_r)
    add_points_in_φ_and_interpolate!(grid, idcs_φ)
    add_points_in_z_and_interpolate!(grid, idcs_z)
end


function extent_in_φ_n_times(grid::CylindricalGrid, n::Int)::CylindricalGrid
    T::Type = eltype(grid.potential)
    l = length(grid.φ)

    g = CylindricalGrid{T}(grid.cyclic, grid.r, zeros(T, l * n), grid.z)
    for i in 1:n
        r = ((i - 1) * l + 1):(i * l)
        g.φ[r] = grid.φ[:] .+ ((i - 1) * grid.cyclic)
        g.potential[:, r, :] = grid.potential[:, :, :]
    end
    return g
end


"""
    extent_2D_grid_to_3D(grid::CylindricalGrid, n::Int)::CylindricalGrid

This function returns a extended grid of a 2-dimensional grid `grid`. The extended grid has
n ticks in φ (symmetrically distributed up to 2π).
"""
function extent_2D_grid_to_3D(grid::CylindricalGrid, n::Int)::CylindricalGrid
    T::Type = eltype(grid.potential)
    l = length(grid.φ)
    if l != 1 error("Grid is not of length 1 -> Not 2D.") end
    cyclic::T = 2π
    φ_step = cyclic / n
    φ::Array{T, 1} = [(i - 1) * φ_step for i in 1:n]
    g = CylindricalGrid{T}(cyclic, grid.r, φ, grid.z)
    for i in 1:n
        r = ((i - 1) * l + 1):(i * l)
        g.potential[:, r, :] = grid.potential[:, :, :]
    end
    return g
end

"""
    extent_2D_grid_to_3D!(grid::CylindricalGrid, n::Int)::Nothing

This function extend a 2-dimensional grid (only one tick in φ) to an 3-dimensional grid with
n ticks in φ (symmetrically distributed up to 2π).
"""
function extent_2D_grid_to_3D!(grid::CylindricalGrid, n::Int)::Nothing
    g = extent_2D_grid_to_3D(grid, n)
    grid.potential = g.potential
    grid.φ = g.φ
    return nothing
end
