struct CylindricalRedBlackGrid{T<:AbstractFloat} <: Grid
    potential::Array{T, 4} # -> last dimension 2 -> red/black

    function CylindricalRedBlackGrid{T}(n_r::Int, n_φ::Int, n_z::Int) where {T<:AbstractFloat}
        potential = zeros(T, n_r, n_φ, n_z, 2)
        return new{T}( potential )
    end
end

function CylindricalRedBlackGrid(grid::CylindricalGrid; only_2d::Bool = false)::CylindricalRedBlackGrid
    odd = 1     
    even = 2 
    T = eltype(grid.potential)
    n_r = length(grid.r) + 2
    n_φ = length(grid.φ) + 2
    nz = length(grid.z)
    n_z = div(nz, 2) + mod(nz, 2) + 2
    cgrb = CylindricalRedBlackGrid{T}( n_r, n_φ, n_z )

    for iz in eachindex(grid.z)
        for iφ in eachindex(grid.φ)
            isum = iz +  iφ
            for ir in eachindex(grid.r)
                if iseven( isum + ir )
                    cgrb.potential[ get_rb_inds(ir, iφ, iz)..., even ] = grid.potential[ ir, iφ, iz ]
                else
                    cgrb.potential[ get_rb_inds(ir, iφ, iz)..., odd ]  = grid.potential[ ir, iφ, iz ]
                end
            end
        end
    end

    # cyclic boundary in φ
    cgrb.potential[:,   1, :, even] = cgrb.potential[:, end - 1, :, even]
    cgrb.potential[:, end, :, even] = cgrb.potential[:,       2, :, even]
    cgrb.potential[:,   1, :, odd]  = cgrb.potential[:, end - 1, :, odd]
    cgrb.potential[:, end, :, odd]  = cgrb.potential[:,       2, :, odd]

    if !only_2d @assert iseven(length(grid.φ)) "Grid is not even in φ => Correct handling of cycling boundary is not implemented yet." end
    return cgrb
end


@inline function get_rb_inds(ir, iφ, iz)::Tuple{Int, Int, Int}
    return ir+1, iφ+1, div(iz, 2) + mod(iz, 2) + 1
end

function get_n_inds(ir, iφ, iz, which::Symbol)::Tuple{Int, Int, Int}
    if which == :even
        if isodd(ir+iφ)
            return ir - 1, iφ - 1, (iz - 1) * 2 - 1
        else
            return ir - 1, iφ - 1, (iz - 1) * 2
        end
    elseif which == :odd
        if iseven(ir+iφ)
            return ir - 1, iφ - 1, (iz - 1) * 2 - 1
        else
            return ir - 1, iφ - 1, (iz - 1) * 2
        end
    else
        error("Symbol 'which' must be ':even' or ':odd'. Given was '$which'.")
    end
end

@inline function get_n_inds(ir, iφ::Int, iz::Int, evenodd::Int)::Tuple{Int, Int, Int}
    if which == :even
        if isodd(ir+iφ)
            return ir - 1, iφ - 1, (iz - 1) * 2 - 1
        else
            return ir - 1, iφ - 1, (iz - 1) * 2
        end
    elseif which == :odd
        if iseven(ir+iφ)
            return ir - 1, iφ - 1, (iz - 1) * 2 - 1
        else
            return ir - 1, iφ - 1, (iz - 1) * 2
        end
    else
        error("Symbol 'which' must be ':even' or ':odd'. Given was '$which'.")
    end
end

function overwrite_potential!(target_grid::CylindricalGrid, source_grid::CylindricalRedBlackGrid)::Nothing
    odd = 1     
    even = 2
    for iz in eachindex(target_grid.z)
        for iφ in eachindex(target_grid.φ)
            sum_zφ = iz + iφ
            for ir in eachindex(target_grid.r)
                target_grid.potential[ir, iφ, iz] = if iseven(ir + sum_zφ)
                    source_grid.potential[get_rb_inds(ir, iφ, iz)..., even]
                else
                    source_grid.potential[get_rb_inds(ir, iφ, iz)..., odd]
                end
            end
        end
    end
    nothing
end

