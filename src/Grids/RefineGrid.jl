
function add_points_in_z_and_interpolate(potential::Array{T, 3}, grid::CylindricalGrid{T}, idcs::Vector{Int})::Tuple{Array{T, 3}, CylindricalGrid{T}} where {T}
    axz::Vector{T} = collect(grid.axes[3])
    n = length(idcs)
    if n == 0 return potential, grid end
    if idcs[end] >= length(axz)
        error("Can not append point in z at the end. Only in between.")
    end

    ranges = [ 1:idcs[1] ]
    for i in (2:n)
        push!(ranges, idcs[i-1]+1:idcs[i])
    end
    push!(ranges, idcs[end]+1:length(axz))

    axz = cat(axz, zeros(T, n), dims=1)

    for (i, range) in enumerate(reverse(ranges))
        axz[range .+ (n+1-i)] = axz[range]
    end
    for (i, range) in enumerate(reverse(ranges[2:end]))
        axz[first(range .+ (n+1-i))-1] = 0.5*(axz[first(range .+ (n+1-i))]+axz[first(range .+ (n+1-i))-2])
    end

    # 1: extend arrays
     potential = cat( potential, zeros(T, size( potential, 1), size( potential, 2), n), dims=3)
    # 2: copy old values to new locations
    for (i, range) in enumerate(reverse(ranges))
        potential[:, :, range .+ (n+1-i)] = potential[:, :, range]
    end
     # 3: interpolate new positions
    for (i, range) in enumerate(reverse(ranges[2:end]))
        potential[:, :, first(range .+ (n+1-i))-1] = 0.5 * ( potential[:, :, first(range .+ (n+1-i))] + potential[:, :, first(range .+ (n+1-i))-2])
    end

    BL::Symbol = typeof(grid.axes[3]).parameters[2]
    BR::Symbol = typeof(grid.axes[3]).parameters[3]
    daz::DiscreteAxis{T, BL, BR} = DiscreteAxis{T, BL, BR}( grid.axes[3].interval, axz )
    grid::CylindricalGrid{T} = CylindricalGrid{T}( (grid.axes[1], grid.axes[2], daz) )

    return potential, grid
end

function add_points_in_r_and_interpolate(potential::Array{T, 3}, grid::CylindricalGrid{T}, idcs::Vector{Int})::Tuple{Array{T, 3}, CylindricalGrid{T}} where {T}
    axr::Vector{T} = collect(grid.axes[1])
    n = length(idcs)
    # idcs = collect(idcs)
    if n == 0 return potential, grid end
    # if !issorted(idcs) idcs = sort(idcs) end
    if idcs[end] >= length(axr)
        error("Can not append point in r at the end. Only in between.")
    end

    ranges = [ 1:idcs[1] ]
    for i in (2:n)
        push!(ranges, idcs[i-1]+1:idcs[i])
    end
    push!(ranges, idcs[end]+1:length(axr))

    axr = cat(axr, zeros(T, n), dims=1)

    for (i, range) in enumerate(reverse(ranges))
        axr[range .+ (n+1-i)] = axr[range]
    end
    for (i, range) in enumerate(reverse(ranges[2:end]))
        axr[first(range .+ (n+1-i))-1] = 0.5*(axr[first(range .+ (n+1-i))]+axr[first(range .+ (n+1-i))-2])
    end

    # 1: extend arrays
    potential = cat( potential, zeros(T, n, size( potential, 2), size( potential, 3)), dims=1)
    # 2: copy old values to new locations
    for (i, range) in enumerate(reverse(ranges))
        potential[range .+ (n+1-i), :, :] = potential[range, :, :]
    end
     # 3: interpolate new positions
    for (i, range) in enumerate(reverse(ranges[2:end]))
        potential[first(range .+ (n+1-i))-1, :, :] = 0.5 * ( potential[first(range .+ (n+1-i)), :, :] + potential[first(range .+ (n+1-i))-2, :, :])
    end

    BL::Symbol = typeof(grid.axes[1]).parameters[2]
    BR::Symbol = typeof(grid.axes[1]).parameters[3]
    da::DiscreteAxis{T, BL, BR} = DiscreteAxis{T, BL, BR}( grid.axes[1].interval, axr )
    grid::CylindricalGrid{T} = CylindricalGrid{T}( (da, grid.axes[2], grid.axes[3]) )

    return potential, grid
end

function add_points_in_φ_and_interpolate(potential::Array{T, 3}, grid::CylindricalGrid{T}, idcs::Vector{Int})::Tuple{Array{T, 3}, CylindricalGrid{T}} where {T}
    axφ::Vector{T} = collect(grid.axes[2])
    n = length(idcs)
    if n == 0 return potential, grid end
    cyclic::T = grid.axes[2].interval.right

    if idcs[end] > length(axφ)
        error("Can not append array in φ at index $(idcs[end]). (Cylcic: $(cyclic)°)")
    end

    append_after_last_index = false
    if idcs[end] == length(axφ)
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
        if idcs[end] < length(axφ)
            push!(ranges, idcs[end]+1:length(axφ))
        end
    end

    if append_after_last_index
        axφ = cat( axφ, zeros(T, n+1), dims=1)
    else
        axφ = cat( axφ, zeros(T, n), dims=1)
    end

    for (i, range) in enumerate(reverse(ranges))
        axφ[range .+ (n+1-i)] = axφ[range]
    end
    for (i, range) in enumerate(reverse(ranges[2:end]))
        axφ[first(range .+ (n+1-i))-1] = 0.5*(axφ[first(range .+ (n+1-i))]+axφ[first(range .+ (n+1-i))-2])
    end

    if append_after_last_index
        axφ[end] = 0.5 * (cyclic + axφ[end-1])
    end

    # 1: extend arrays
    if append_after_last_index
        potential = cat( potential, zeros(T, size( potential, 1), n + 1, size( potential, 3)), dims=2 )
    else
        potential = cat( potential, zeros(T, size( potential, 1), n, size( potential, 3)), dims=2)
    end
    # 2: copy old values to new locations
    for (i, range) in enumerate(reverse(ranges))
        potential[:, range .+ (n+1-i), :] = potential[:, range, :]
    end
     # 3: interpolate new positions
    for (i, range) in enumerate(reverse(ranges[2:end]))
        potential[: ,first(range .+ (n+1-i))-1, :] = 0.5 * ( potential[:, first(range .+ (n+1-i)), :] + potential[:, first(range .+ (n+1-i))-2, :])
    end
    if append_after_last_index
        potential[: ,end, :] = 0.5 * ( potential[:, end-1, :] + potential[:, 1, :])
    end

    BL::Symbol = typeof(grid.axes[2]).parameters[2]
    BR::Symbol = typeof(grid.axes[2]).parameters[3]
    da::DiscreteAxis{T, BL, BR} = DiscreteAxis{T, BL, BR}( grid.axes[2].interval, axφ )
    grid::CylindricalGrid{T} = CylindricalGrid{T}( (grid.axes[1], da, grid.axes[3]) )

    return potential, grid
end


function add_points_and_interpolate(potential::Array{T, 3}, grid::CylindricalGrid{T}, idcs_r::Vector{Int}, idcs_φ::Vector{Int}, idcs_z::Vector{Int})::Tuple{Array{T, 3}, CylindricalGrid{T}} where {T}
    potential::Array{T, 3}, grid::CylindricalGrid{T} = add_points_in_r_and_interpolate(potential, grid, idcs_r)
    potential, grid = add_points_in_φ_and_interpolate(potential, grid, idcs_φ)
    if isodd(length(grid.axes[2])) && length(length(grid.axes[2])) != 1
        potential, grid = add_points_in_φ_and_interpolate(potential, grid, [1])
    end    
    potential, grid = add_points_in_z_and_interpolate(potential, grid, idcs_z)
    potential, grid
end

function add_points_and_interpolate(potential::Array{T, 3}, grid::CartesianGrid3D{T}, idcs_x::Vector{Int}, idcs_y::Vector{Int}, idcs_z::Vector{Int})::Tuple{Array{T, 3}, CartesianGrid3D{T}} where {T}
    potential::Array{T, 3}, grid::CartesianGrid3D{T} = add_points_in_x_and_interpolate(potential, grid, idcs_x)
    potential, grid = add_points_in_y_and_interpolate(potential, grid, idcs_y)
    potential, grid = add_points_in_z_and_interpolate(potential, grid, idcs_z)
    potential, grid
end


function add_points_in_z_and_interpolate(potential::Array{T, 3}, grid::CartesianGrid3D{T}, idcs::Vector{Int})::Tuple{Array{T, 3}, CartesianGrid3D{T}} where {T}
    ax::Vector{T} = collect(grid.axes[3])
    n = length(idcs)
    if n == 0 return potential, grid end
    if idcs[end] >= length(ax)
        error("Can not append point in z at the end. Only in between.")
    end

    ranges = [ 1:idcs[1] ]
    for i in (2:n)
        push!(ranges, idcs[i-1]+1:idcs[i])
    end
    push!(ranges, idcs[end]+1:length(ax))

    ax = cat(ax, zeros(T, n), dims = 1)

    for (i, range) in enumerate(reverse(ranges))
        ax[range .+ (n + 1 - i)] = ax[range]
    end
    for (i, range) in enumerate(reverse(ranges[2:end]))
        ax[first(range .+ (n + 1 - i)) - 1] = 0.5 * (ax[first(range .+ (n + 1 - i))] + ax[first(range .+ (n + 1 - i)) - 2])
    end

    # 1: extend arrays
    potential = cat( potential, zeros(T, size( potential, 1), size( potential, 2), n), dims = 3)
    # 2: copy old values to new locations
    for (i, range) in enumerate(reverse(ranges))
        potential[:, :, range .+ (n+1-i)] = potential[:, :, range]
    end
     # 3: interpolate new positions
    for (i, range) in enumerate(reverse(ranges[2:end]))
        potential[:, :, first(range .+ (n+1-i))-1] = 0.5 * ( potential[:, :, first(range .+ (n+1-i))] + potential[:, :, first(range .+ (n+1-i))-2])
    end

    BL::Symbol = typeof(grid.axes[3]).parameters[2]
    BR::Symbol = typeof(grid.axes[3]).parameters[3]
    da::DiscreteAxis{T, BL, BR} = DiscreteAxis{T, BL, BR}( grid.axes[3].interval, ax )
    grid::CartesianGrid3D{T} = CartesianGrid3D{T}( (grid.axes[1], grid.axes[2], da) )

    return potential, grid
end

function add_points_in_x_and_interpolate(potential::Array{T, 3}, grid::CartesianGrid3D{T}, idcs::Vector{Int})::Tuple{Array{T, 3}, CartesianGrid3D{T}} where {T}
    ax::Vector{T} = collect(grid.axes[1])
    n = length(idcs)
    if n == 0 return potential, grid end
    if idcs[end] >= length(ax)
        error("Can not append point in z at the end. Only in between.")
    end

    ranges = [ 1:idcs[1] ]
    for i in (2:n)
        push!(ranges, idcs[i-1]+1:idcs[i])
    end
    push!(ranges, idcs[end]+1:length(ax))

    ax = cat(ax, zeros(T, n), dims = 1)

    for (i, range) in enumerate(reverse(ranges))
        ax[range .+ (n + 1 - i)] = ax[range]
    end
    for (i, range) in enumerate(reverse(ranges[2:end]))
        ax[first(range .+ (n + 1 - i)) - 1] = 0.5 * (ax[first(range .+ (n + 1 - i))] + ax[first(range .+ (n + 1 - i)) - 2])
    end

    # 1: extend arrays
    potential = cat( potential, zeros(T, n, size( potential, 2), size( potential, 3)), dims = 1)
    # 2: copy old values to new locations
    for (i, range) in enumerate(reverse(ranges))
        potential[range .+ (n+1-i), :, :] = potential[range, :, :]
    end
     # 3: interpolate new positions
    for (i, range) in enumerate(reverse(ranges[2:end]))
        potential[first(range .+ (n+1-i))-1, :, :] = 0.5 * ( potential[first(range .+ (n+1-i)), :, :] + potential[first(range .+ (n+1-i))-2, :, :])
    end

    BL::Symbol = typeof(grid.axes[1]).parameters[2]
    BR::Symbol = typeof(grid.axes[1]).parameters[3]
    da::DiscreteAxis{T, BL, BR} = DiscreteAxis{T, BL, BR}( grid.axes[1].interval, ax )
    grid::CartesianGrid3D{T} = CartesianGrid3D{T}( (da, grid.axes[2], grid.axes[3]) )

    return potential, grid
end

function add_points_in_y_and_interpolate(potential::Array{T, 3}, grid::CartesianGrid3D{T}, idcs::Vector{Int})::Tuple{Array{T, 3}, CartesianGrid3D{T}} where {T}
    ax::Vector{T} = collect(grid.axes[2])
    n = length(idcs)
    if n == 0 return potential, grid end
    if idcs[end] >= length(ax)
        error("Can not append point in z at the end. Only in between.")
    end

    ranges = [ 1:idcs[1] ]
    for i in (2:n)
        push!(ranges, idcs[i-1]+1:idcs[i])
    end
    push!(ranges, idcs[end]+1:length(ax))

    ax = cat(ax, zeros(T, n), dims = 1)

    for (i, range) in enumerate(reverse(ranges))
        ax[range .+ (n + 1 - i)] = ax[range]
    end
    for (i, range) in enumerate(reverse(ranges[2:end]))
        ax[first(range .+ (n + 1 - i)) - 1] = 0.5 * (ax[first(range .+ (n + 1 - i))] + ax[first(range .+ (n + 1 - i)) - 2])
    end

    # 1: extend arrays
    potential = cat( potential, zeros(T, size( potential, 1), n, size( potential, 3)), dims = 2)
    # 2: copy old values to new locations
    for (i, range) in enumerate(reverse(ranges))
        potential[:, range .+ (n+1-i), :] = potential[:, range, :]
    end
     # 3: interpolate new positions
    for (i, range) in enumerate(reverse(ranges[2:end]))
        potential[:, first(range .+ (n+1-i))-1, :] = 0.5 * ( potential[:, first(range .+ (n+1-i)), :] + potential[:, first(range .+ (n+1-i))-2, :])
    end

    BL::Symbol = typeof(grid.axes[2]).parameters[2]
    BR::Symbol = typeof(grid.axes[2]).parameters[3]
    da::DiscreteAxis{T, BL, BR} = DiscreteAxis{T, BL, BR}( grid.axes[2].interval, ax )
    grid::CartesianGrid3D{T} = CartesianGrid3D{T}( (grid.axes[1], da,  grid.axes[3]) )

    return potential, grid
end