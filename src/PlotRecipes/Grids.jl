
@recipe function f(grid::Grid{T, N, S}) where {T, N, S}    
    #only plot grid point densities for axes with more than one tick
    iaxs = findall(ax -> length(ax) > 1, grid.axes)
    layout --> length(iaxs) 
    seriestype --> :stephist
    legend --> false

    for (i,iax) in enumerate(iaxs)
        @series begin
            subplot := i
            bins --> div(length(grid[iax]), 2)
            xguide --> "Grid point density - Axis $(iax)"
            grid[iax]
        end
    end
end
