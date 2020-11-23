@recipe function f(det::SolidStateDetector{T}; SSD_style = :wireframe, n = 30, φ = missing, seriescolor = missing, label = missing, alpha_factor = 1) where {T}
    if !(SSD_style in [:wireframe, :samplesurface])
        @warn "Chose SSD_style from [:wireframe, :samplesurface]. Defaulting to :wireframe"
        SSD_style = :wireframe
    end
    clabel = (ismissing(label) ? map(c -> (c.name == "" ? c.id : c.name), det.contacts) : label)
    if !(typeof(clabel) <: AbstractArray) clabel = [clabel] end
    ccolor = (ismissing(seriescolor) ? map(c -> c.id, det.contacts) : seriescolor)
    if !(typeof(ccolor) <: AbstractArray) ccolor = [ccolor] end

    if ismissing(φ)
        xguide --> "x / m"
        yguide --> "y / m"
        zguide --> "z / m"

        for (cn,contact) in enumerate(det.contacts)
            linewidth --> 2
            seriescolor := ccolor[(cn-1)%length(ccolor)+1]
            @series begin
                label := ""
                n --> n
                SSD_style --> SSD_style
                detector --> det
                alpha_factor --> alpha_factor
                contact
            end
            @series begin
                label := clabel[(cn-1)%length(clabel)+1]
                seriescolor := ccolor[(cn-1)%length(ccolor)+1]
                []
            end
        end
    else
        @info "2D plotting of the detector not yet implemented."
    end
end


@recipe function f(contact::Contact{T}; SSD_style = :wireframe, n = 30, seriescolor = missing, detector = missing, alpha_factor = 1) where {T}
    ccolor = (ismissing(seriescolor) ? contact.id : seriescolor)
    @series begin
        seriescolor := ccolor
        label --> (contact.name == "" ? contact.id : contact.name)
        []
    end
    for (cn,c) in enumerate(contact.geometry_positive)
        @series begin
            seriescolor := ccolor
            label := ""
            n --> n
            SSD_style --> SSD_style
            detector --> detector
            contact --> contact
            alpha_factor --> alpha_factor
            c
        end
    end
end

#=
@recipe function f(geometry::AbstractGeometry)
    pos, neg = get_decomposed_volumes(geometry)
    for g in pos
        @series begin
            g
        end
    end
end
@recipe function f(object::AbstractObject)
    @series begin
        object.geometry
    end
end

## 2D from Felix


function line_2d(r1,r2,z1,z2)
   x1::Vector = [r1,z1]
   x2::Vector = [r2,z2]
   return [x1[1],x2[1]],[x1[2],x2[2]]
end


function polar_circle_2d(r::T, phiStart::T, phiStop::T;nSteps=360) where T <: SSDFloat
    if r == 0; return [],[]; end
    phiRange = collect(phiStart:(phiStop-phiStart)/nSteps:phiStop)
    rRange = [r for a in 1:length(phiRange)]
    return phiRange, rRange
end




@recipe function f(d::SolidStateDetector{T}, dim::Symbol; φ = missing, z = missing) where{T <: SSDFloat}
    if ismissing(z) && ismissing(φ)
        print("Please specify φ or z.")
    elseif ismissing(z)
        for c in d.contacts
            for g in c.geometry_positive
                @series begin
                    if d.name == "Public Inverted Coax"
                        if typeof(c) == Contact{T,:N}; color --> :orange
                        elseif typeof(c) == Contact{T,:P}; color --> :blue; end
                    else color --> :black
                    end
                    lw --> 2
                    label = ""
                    g, :φ, T(deg2rad(mod(φ,2π)))
                end
            end
        end
    elseif ismissing(φ)
        proj --> :polar
        for c in d.contacts
            for g in c.geometry_positive
                @series begin
                    if d.name == "Public Inverted Coax"
                        if typeof(c) == Contact{T,:N}; color --> :orange
                        elseif typeof(c) == Contact{T,:P}; color --> :blue; end
                    else color --> :black
                    end
                    lw --> 2
                    label = ""
                    g, :z, T(z)
                end
            end
        end
    else
        print("Please specify only one of either φ or z.")
    end
end


@recipe function f(Vol::Tube{T}, dim::Symbol, parameter::T) where{T <: SSDFloat}
    if dim == :φ
        if parameter in Vol.φ_interval
           rStart = Vol.r_interval.left
           rStop = Vol.r_interval.right
           zStart = Vol.z_interval.left
           zStop = Vol.z_interval.right

           @series begin
                   label --> ""
                   line_2d(rStop,rStop,zStart,zStop)
               end

            @series begin
                   label --> ""
                   line_2d(rStart,rStop,zStart,zStart)
            end

           if rStart != rStop
                @series begin
                   label --> ""
                   line_2d(rStart,rStop,zStop,zStop)
               end
            end

            if zStart != zStop && rStart != 0
                @series begin
                   label --> ""
                   line_2d(rStart,rStart,zStart,zStop)
               end
            end
        end

    elseif dim == :z
        if parameter in Vol.z_interval
            rStart = Vol.r_interval.left
            rStop = Vol.r_interval.right
            phiStart = Vol.φ_interval.left
            phiStop = Vol.φ_interval.right
            zStart = Vol.z_interval.left
            zStop = Vol.z_interval.right

            if rStart != 0
                @series begin
                    label --> ""
                    polar_circle_2d(rStart,phiStart,phiStop)
                end
            end
            @series begin
                label --> ""
                polar_circle_2d(rStop,phiStart,phiStop)
            end

            if !isapprox(mod(phiStart,2π),mod(phiStop,2π),atol = 1e-5)
                @series begin
                    label --> ""
                    line_2d(phiStart,phiStart,rStart,rStop)
                end
                @series begin
                    label --> ""
                    line_2d(phiStop,phiStop,rStart,rStop)
                end
            end
        end
    end
end



# @recipe function f(Vol::ConeMantle{T}, dim::Symbol, parameter::T) where{T <: SSDFloat}
#     newVol = Vol.cone
#     @series begin
#         newVol, dim, parameter
#     end
# end



@recipe function f(Vol::Cone{T}, dim::Symbol, parameter::T) where{T <: SSDFloat}
    if dim == :φ
        if parameter in Vol.φ_interval
            rStart = Vol.r_interval.left
            rStop = Vol.r_interval.right
            zStart = Vol.z_interval.left
            zStop = Vol.z_interval.right
            orientation = Vol.orientation

            if orientation in [:top_right, :bottom_left]
                @series begin
                    label --> ""
                   line_2d(rStop,rStart,zStart,zStop)
                end
            else
                @series begin
                    label --> ""
                   line_2d(rStart,rStop,zStart,zStop)
                end
            end
        end

    elseif dim == :z
        if parameter in Vol.z_interval
            rStart = Vol.r_interval.left
            rStop = Vol.r_interval.right
            phiStart = Vol.φ_interval.left
            phiStop = Vol.φ_interval.right
            zStart = Vol.z_interval.left
            zStop = Vol.z_interval.right
            orientation = Vol.orientation

            if orientation in [:top_right, :bottom_left]
                @series begin
                    label --> ""
                    r = (rStop-rStart)*(parameter-zStop)/(zStart-zStop) + rStart
                    polar_circle_2d(r,phiStart,phiStop)
                end
            else
                @series begin
                    label --> ""
                    r = (rStop-rStart)*(parameter-zStart)/(zStop-zStart) + rStart
                    polar_circle_2d(r,phiStart,phiStop)
                end
            end
        end

    end
end
=#
