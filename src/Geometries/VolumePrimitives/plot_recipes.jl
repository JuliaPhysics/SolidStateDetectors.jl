@recipe function f(Vol::SSD.Tube{T}, coloring = missing) where T <: SSDFloat
    rStart = Vol.r_interval.left
    rStop = Vol.r_interval.right
    φStart = Vol.φ_interval.left
    φStop = Vol.φ_interval.right
    zStart = Vol.z_interval.left
    zStop = Vol.z_interval.right

    ismissing(Vol.translate) ? translate = [0.0,0.0,0.0] : translate = Vol.translate
    ismissing(coloring) ? c --> :green : c --> coloring
    @series begin
        label --> "Tube"
        partialcircle_3d(rStop,φStart,φStop,[0,0,zStart] .+ translate)
    end
    @series begin
        label:= ""
        partialcircle_3d(rStop,φStart,φStop,[0,0,zStop] .+ translate)
    end
    @series begin
        label:= ""
        partialcircle_3d(rStart,φStart,φStop,[0,0,zStart] .+ translate)
    end
    @series begin
        label:= ""
        partialcircle_3d(rStart,φStart,φStop,[0,0,zStop] .+ translate)
    end
    ## Vertical Lines

    if !iszero(rStart)
        @series begin
            label:= ""
            line_3d(rStart, rStart, φStart, φStart, zStart, zStop, translate = translate)
        end
    end
    if !iszero(rStart)
        @series begin
            label:= ""
            line_3d(rStart, rStart, φStop, φStop, zStart, zStop, translate = translate)
        end
    end
    @series begin
        label:= ""
        line_3d(rStop, rStop, φStart, φStart, zStart, zStop, translate = translate)
    end
    @series begin
        label:= ""
        line_3d(rStop, rStop, φStop, φStop, zStart, zStop, translate = translate)
    end

    ##Horizontal Lines

    if !isapprox((φStop - φStart)%2π , 0.0,atol=0.00001)
        @series begin
            label:= ""
            line_3d(rStart, rStop, φStart, φStart, zStart, zStart, translate = translate)
        end
        @series begin
            label:= ""
            line_3d(rStart, rStop, φStop, φStop, zStart, zStart, translate = translate)
        end
        @series begin
            label:= ""
            line_3d(rStart, rStop, φStart, φStart, zStop, zStop, translate = translate)
        end
        @series begin
            label:= ""
            line_3d(rStart, rStop, φStop, φStop, zStop, zStop, translate = translate)
        end
    end
end


@recipe function f(vol::SSD.Cone{T},n_aux_lines = 0, coloring = missing) where T
    ismissing(vol.translate) ? translate = [0.0,0.0,0.0] : translate = vol.translate
    ismissing(coloring) ? c --> :orange : c --> coloring
    @series begin
        label --> "Cone"
        partialcircle_3d(vol.rStop1,vol.φStart,vol.φStop,[0,0,vol.zStart]+translate)
    end
    @series begin
        label:= ""
        partialcircle_3d(vol.rStop2,vol.φStart,vol.φStop,[0,0,vol.zStop]+translate)
    end
    @series begin
        label:= ""
        partialcircle_3d(vol.rStart1,vol.φStart,vol.φStop,[0,0,vol.zStart]+translate)
    end
    @series begin
        label:= ""
        c-->coloring
        partialcircle_3d(vol.rStart2,vol.φStart,vol.φStop,[0,0,vol.zStop]+translate)
    end

    ## Vertical Lines


    @series begin
        label:= ""
        line_3d(vol.rStart1,vol.rStart2,vol.φStart,vol.φStart,vol.zStart,vol.zStop, translate = translate)
    end

    @series begin
        label:= ""
        line_3d(vol.rStart1,vol.rStart2,vol.φStop,vol.φStop,vol.zStart,vol.zStop, translate = translate)
    end


    @series begin
        label:= ""
        line_3d(vol.rStop1,vol.rStop2,vol.φStart,vol.φStart,vol.zStart,vol.zStop, translate = translate)
    end
    @series begin
        label:= ""
        line_3d(vol.rStop1,vol.rStop2,vol.φStop,vol.φStop,vol.zStart,vol.zStop, translate = translate)
    end

    ##Horizontal Lines

    if !isapprox((vol.φStop - vol.φStart)%2π , 0.0,atol=0.00001)
        @series begin
            label:= ""
            line_3d(vol.rStart1,vol.rStop1,vol.φStart,vol.φStart,vol.zStart,vol.zStart, translate = translate)
        end
        @series begin
            label:= ""
            line_3d(vol.rStart1,vol.rStop1,vol.φStop,vol.φStop,vol.zStart,vol.zStart, translate = translate)
        end
        @series begin
            label:= ""
            line_3d(vol.rStart2,vol.rStop2,vol.φStart,vol.φStart,vol.zStop,vol.zStop, translate = translate)
        end
        @series begin
            label:= ""
            line_3d(vol.rStart2,vol.rStop2,vol.φStop,vol.φStop,vol.zStop,vol.zStop, translate = translate)
        end
    end

end
function partialcircle_3d(radius,phiStart,phiStop,Translate::AbstractVector;nSteps=400)
    phirange = mylinspace(phiStart,phiStop,nSteps)

    x::Vector{AbstractFloat}=map(x->radius*cos.(x) , phirange)
    y::Vector{AbstractFloat}=map(x->radius*sin.(x) , phirange)
    z::Vector{AbstractFloat}=map(x->Translate[3], phirange)
    x.+=Translate[1]
    y.+=Translate[2]

    return x,y,z
end
