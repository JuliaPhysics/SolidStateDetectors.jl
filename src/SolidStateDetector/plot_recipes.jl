@recipe function f(::Val{:InvertedCoax}, dim::Symbol;
                    r=missing,
                    φ=missing,
                    z=missing,
                    half=true)
    T = typeof(d.crystal_radius)

    cross_section::Symbol, value::T = if ismissing(φ) && ismissing(r) && ismissing(z)
        :φ, 0
    elseif !ismissing(φ) && ismissing(r) && ismissing(z)
        :φ, φ
    elseif ismissing(φ) && !ismissing(r) && ismissing(z)
        :r, r
    elseif ismissing(φ) && ismissing(r) && !ismissing(z)
        :z, z
    else
        error(ArgumentError, ": Only one of the keywords `r, φ, z` is allowed.")
    end

    if cross_section == :φ
        points_to_connect::AbstractVector = [SVector{2,T}(0.0,0.0)]
        add_point!(points_to_connect, d.groove_rInner,T(0.0))
        add_point!(points_to_connect, d.groove_rInner,d.groove_depth)
        add_point!(points_to_connect, d.groove_rInner+d.groove_width,d.groove_depth)
        add_point!(points_to_connect, d.groove_rInner+d.groove_width,T(0.0))
        add_point!(points_to_connect, d.crystal_radius,T(0.0))
        add_point!(points_to_connect, d.crystal_radius,d.crystal_length-d.taper_outer_length)
        add_point!(points_to_connect, d.taper_outer_rInner,d.crystal_length)
        add_point!(points_to_connect, d.borehole_radius,d.crystal_length)
        add_point!(points_to_connect, d.borehole_radius,d.crystal_length-d.borehole_length)
        add_point!(points_to_connect, T(0.0),d.crystal_length-d.borehole_length)
        if half == false
            for ipoint in points_to_connect[end-1:-1:1 ]
                add_point!(points_to_connect,-1*ipoint[1],ipoint[2])
            end
        end
        @series begin
            lw --> 1
            label --> ""
            linestyle --> :dash
            color --> :red
            points_to_connect
        end
    end
end

@recipe function f(::Val{:BEGe}, dim::Symbol;
                    r=missing,
                    φ=missing,
                    z=missing,
                    half=true)
    T = typeof(d.crystal_radius)

    cross_section::Symbol, value::T = if ismissing(φ) && ismissing(r) && ismissing(z)
        :φ, 0
    elseif !ismissing(φ) && ismissing(r) && ismissing(z)
        :φ, φ
    elseif ismissing(φ) && !ismissing(r) && ismissing(z)
        :r, r
    elseif ismissing(φ) && ismissing(r) && !ismissing(z)
        :z, z
    else
        error(ArgumentError, ": Only one of the keywords `r, φ, z` is allowed.")
    end

    if cross_section == :φ
        points_to_connect::AbstractVector = [SVector{2,T}(0.0,0.0)]
        add_point!(points_to_connect, d.groove_rInner,T(0.0))
        if d.groove_endplate == "bot"
            add_point!(points_to_connect, d.groove_rInner,T(0.0))
            add_point!(points_to_connect, d.groove_rInner,T(0.0)+d.groove_depth)
            add_point!(points_to_connect, d.groove_rInner+d.groove_width,T(0.0)+d.groove_depth)
            add_point!(points_to_connect, d.groove_rInner+d.groove_width,T(0.0))
        end
        add_point!(points_to_connect, d.taper_bot_rInner,T(0.0))
        add_point!(points_to_connect, d.crystal_radius,d.taper_bot_length)
        add_point!(points_to_connect, d.crystal_radius,d.crystal_length-d.taper_top_length)
        add_point!(points_to_connect, d.taper_top_rInner,d.crystal_length)
        if d.groove_endplate == "top"
            add_point!(points_to_connect, d.groove_rInner+d.groove_width,d.crystal_length)
            add_point!(points_to_connect, d.groove_rInner+d.groove_width,d.crystal_length-d.groove_depth)
            add_point!(points_to_connect, d.groove_rInner,d.crystal_length-d.groove_depth)
            add_point!(points_to_connect, d.groove_rInner,d.crystal_length)
        end
        add_point!(points_to_connect, T(0.0),d.crystal_length)
        if half == false
            for ipoint in points_to_connect[end-1:-1:1 ]
                add_point!(points_to_connect,-1*ipoint[1],ipoint[2])
            end
        end
        @series begin
            lw --> 1
            label --> ""
            linestyle --> :dash
            color --> :red
            points_to_connect
        end
    end
end

function add_point!(points_to_connect::AbstractVector{<:SVector{2,T}},x::T,y::T) where T <: Real
    push!(points_to_connect,SVector{2,T}(x,y))
end
@recipe function f(::Val{:Coax})
    size:=(800,800)
    linewidth:=2
    linecolor:=:black
    label-->""
    # aspect_ratio := equal
    f=1/c.geometry_unit_factor
    @series begin
        partialcircle_3d(f*c.crystal_radius,0,2π,[0,0,0])
    end
    @series begin
        partialcircle_3d(f*c.crystal_radius,0,2π,[0,0,f*c.crystal_length])
    end
    @series begin
        partialcircle_3d(f*c.taper_inner_top_rOuter,0,2π,[0,0,f*c.crystal_length])
    end
    linestyle:=:dash
    @series begin
        partialcircle_3d(f*c.borehole_radius,0,2π,[0,0,f*c.crystal_length-f*c.taper_inner_top_length])
    end
    @series begin
        partialcircle_3d(f*c.taper_inner_bot_rOuter,0,2π,[0,0,0])
    end
    @series begin
        partialcircle_3d(f*c.borehole_radius,0,2π,[0,0,f*c.taper_inner_bot_length])
    end

    for i_rep_seg in 2:size(c.segmentation_r_ranges,1)
        @series begin
            partialcircle_3d(f*c.crystal_radius,c.segmentation_phi_ranges[i_rep_seg][1],c.segmentation_phi_ranges[i_rep_seg][2],[0,0,c.segmentation_z_ranges[i_rep_seg][1]*f])
        end
        @series begin
            partialcircle_3d(f*c.crystal_radius,c.segmentation_phi_ranges[i_rep_seg][1],c.segmentation_phi_ranges[i_rep_seg][2],[0,0,c.segmentation_z_ranges[i_rep_seg][2]*f])
        end
        @series begin
            x::Vector{AbstractFloat}=[f*c.crystal_radius*cos.(c.segmentation_phi_ranges[i_rep_seg][1]),f*c.crystal_radius*cos.(c.segmentation_phi_ranges[i_rep_seg][1])]
            y::Vector{AbstractFloat}=[f*c.crystal_radius*sin.(c.segmentation_phi_ranges[i_rep_seg][1]),f*c.crystal_radius*sin.(c.segmentation_phi_ranges[i_rep_seg][1])]
            z::Vector{AbstractFloat}=[c.segmentation_z_ranges[i_rep_seg][1]*f,c.segmentation_z_ranges[i_rep_seg][2]*f]
            x,y,z
        end
        @series begin
            x=[f*c.crystal_radius*cos.(c.segmentation_phi_ranges[i_rep_seg][2]),f*c.crystal_radius*cos.(c.segmentation_phi_ranges[i_rep_seg][2])]
            y=[f*c.crystal_radius*sin.(c.segmentation_phi_ranges[i_rep_seg][2]),f*c.crystal_radius*sin.(c.segmentation_phi_ranges[i_rep_seg][2])]
            z=[c.segmentation_z_ranges[i_rep_seg][1]*f,c.segmentation_z_ranges[i_rep_seg][2]*f]
            x,y,z
        end
    end
end

struct outer_taper{T <: SSDFloat}
    rStart::T
    rStop::T
    φStart::T
    φStop::T
    zStart::T
    zStop::T
    orientation::String
    function outer_taper(rStart::T, rStop::T, φStart::T, φStop::T, zStart::T, zStop::T, orientation) where T <: SSDFloat
        return new{T}(rStart, rStop, φStart, φStop, zStart, zStop, orientation)
    end
end

function outer_taper(rStart, rStop, φStart, φStop, zStart, zStop, orientation)
    return outer_taper{typeof(rStart)}(rStart, rStop, φStart, φStop, zStart, zStop, orientation)
end

@recipe function f(::Val{:InvertedCoax}; coloring=[], labeling=[])
    if b.name =="ExampleInvertedCoax"
        coloring = [:blue, :orange, :orange, :orange, :orange, :orange, :orange, :orange]
        labeling = ["Core", "Mantle", "", "", "", "", "", ""]
    end
    # aspect_ratio := equal
    f=1/b.geometry_unit_factor
    for i_part in 1:size(b.segmentation_r_ranges,1)
        rStart = f*b.segmentation_r_ranges[i_part][1]
        rStop = f*b.segmentation_r_ranges[i_part][2]
        ΘStart = b.segmentation_phi_ranges[i_part][1]
        ΘStop = b.segmentation_phi_ranges[i_part][2]
        zStart = f*b.segmentation_z_ranges[i_part][1]
        zStop = f*b.segmentation_z_ranges[i_part][2]
        if !isempty(coloring)
            color := coloring[i_part]
        else
            color := :black
        end
        if b.segmentation_types[i_part]=="Tubs"
            if !isempty(labeling)
                label := labeling[i_part]
            else
                label := ""
            end
            if b.borehole_modulation == true && i_part in [b.borehole_segment_idx,b.borehole_bot_segment_idx]
                @series begin
                    phirange = [i for i in 0:0.05:2π+0.5]
                    xs = [(b.borehole_ModulationFunction(phirange[i])+b.borehole_radius) * cos(phirange[i]) for i in eachindex(phirange)] .* f
                    ys = [(b.borehole_ModulationFunction(phirange[i])+b.borehole_radius)* sin(phirange[i]) for i in eachindex(phirange)] .* f
                    zs = [zStart for i in eachindex(phirange)]
                    xs, ys, zs
                end
            elseif b.borehole_modulation == true && i_part == b.borehole_top_segment_idx
                @series begin
                    phirange = [i for i in 0:0.05:2π+0.5]
                    xs = [(b.borehole_ModulationFunction(phirange[i])+b.borehole_radius) * cos(phirange[i]) for i in eachindex(phirange)] .* f
                    ys = [(b.borehole_ModulationFunction(phirange[i])+b.borehole_radius)* sin(phirange[i]) for i in eachindex(phirange)] .* f
                    zs = [zStart for i in eachindex(phirange)]
                    xs, ys, zs
                end
                @series begin
                    phirange = [i for i in 0:0.05:2π+0.5]
                    xs = [rStop * cos(phirange[i]) for i in eachindex(phirange)]
                    ys = [rStop * sin(phirange[i]) for i in eachindex(phirange)]
                    zs = [zStart for i in eachindex(phirange)]
                    xs, ys, zs
                end
            else
                @series begin
                    Tubs("",1,1.,"",0.0,rStart,rStop,ΘStart,ΘStop,zStart,zStop)
                end
            end
        else
            @series begin
                if !isempty(labeling)
                    label := labeling[i_part]
                else
                    label := ""
                end
                # Taper("",1,1.,"",0.0,rStart,rStop,ΘStart,ΘStop,zStart,zStop,"bl")
                outer_taper(rStart,rStop,ΘStart,ΘStop,zStart,zStop,"c//")
            end
        end
    end
end

@recipe function f(::Val{:BEGe}; coloring=[], labeling=[])
    if b.name =="ExampleSegmentedBEGe"
        coloring = [:blue, :red, :purple, :orange, :purple, :grey, :purple,
                            :red, :purple, :orange, :purple, :grey, :purple,
                            :red, :purple, :orange, :purple, :grey, :purple, :purple]
        labeling = ["Core", "Seg. 1", "", "Seg. 2", "", "Seg.3", "Seg. 4",
                            "","","","","","",
                            "","","","","","",""]
        xlims --> (-45,45)
        ylims --> (-45,45)
        zlims --> (-25,65)
    end
    # aspect_ratio := equal
    f=1/b.geometry_unit_factor
    for i_part in 1:size(b.segmentation_r_ranges,1)
        rStart = f*b.segmentation_r_ranges[i_part][1]
        rStop = f*b.segmentation_r_ranges[i_part][2]
        φStart = b.segmentation_phi_ranges[i_part][1]
        φStop = b.segmentation_phi_ranges[i_part][2]
        zStart = f*b.segmentation_z_ranges[i_part][1]
        zStop = f*b.segmentation_z_ranges[i_part][2]
        if !isempty(coloring)
            color := coloring[i_part]
        else
            color := :black
        end
        if b.segmentation_types[i_part]=="Tubs"
            @series begin
                if !isempty(labeling)
                    label := labeling[i_part]
                else
                    label := ""
                end
                Tubs("",1,1.,"",0.0,rStart,rStop,φStart,φStop,zStart,zStop)
            end
        else
            @series begin
                if !isempty(labeling)
                    label := labeling[i_part]
                else
                    label := ""
                end
                # Taper("",1,1.,"",0.0,rStart,rStop,φStart,φStop,zStart,zStop,"bl")
                outer_taper(rStart,rStop,φStart,φStop,zStart,zStop,"c//")
            end
        end
    end
end

function mylinspace(Start,Stop,nSteps)
    return collect(Start:(Stop-Start)/nSteps:Stop)
end



function partialcircle(radius,phiStart,phiStop,Translate::Vector=[0.0,0.0,0.0];nSteps=400)
    phirange = mylinspace(phiStart,phiStop,nSteps)
    x::Vector{AbstractFloat}=map(x->radius*cos.(x),phirange)
    y::Vector{AbstractFloat}=map(x->radius*sin.(x),phirange)
    x.+=Translate[1]
    y.+=Translate[2]
    return x,y
end

function line_3d(r1::T, r2::T, phi1::T, phi2::T, z1::T, z2::T; translate::AbstractVector = SVector{3,T}(0.0,0.0,0.0)) where T <:Real
    x1::SVector{3,T} = [r1*cos(phi1),r1*sin(phi1),z1]+translate
    x2::SVector{3,T} = [r2*cos(phi2),r2*sin(phi2),z2]+translate
    return [x1[1],x2[1]] , [x1[2],x2[2]], [x1[3],x2[3]]
end

@recipe function f(contact::AbstractContact{T}) where T
    c-->:orange
    for (i,g) in enumerate(vcat(contact.geometry_positive, contact.geometry_negative))
        @series begin
            i==1 ? label --> "$(contact.id)" : label := ""
            g
        end
    end
end

@recipe function f(c::SolidStateDetector{T}; coloring = [], labeling = [], gfs = 13, lfs = 16, tfs = 13) where T
    legendfontsize := lfs
    # aspect_ratio := 1
    tickfontsize := tfs
    guidefontsize := gfs
    lw --> 2
    # xlabel --> "\n x / m"
    # ylabel --> "\n y / m"
    # zlabel --> "\n z / m"
    camera --> (15,25)
    if ismissing(coloring) coloring = 1:length(c.contacts) end
    if ismissing(labeling) labeling = ["ChnID=$(contact.id)" for contact in c.contacts] end
    labeling = labeling
    # if c.name == "Public Segmented BEGe"
    #     @series begin
    #         # zlims --> (-0.023,0.063)
    #         # aspect_ratio --> 0.5
    #         c, Val(:segBEGe)
    #     end
    # elseif c.name == "Public Inverted Coax"
    #     @series begin
    #         c, Val(:ivc)
    #     end
    # elseif c.name == "Public Coax"
    #     @series begin
    #         c, Val(:coax)
    #     end
    # else
    begin
        for (ic,contact) in enumerate(c.contacts)
            @series begin
                guidefontsize --> 22
                !isempty(coloring) ? c --> coloring[ic] : nothing
                !isempty(labeling) ? label --> labeling[ic] : nothing
                contact
            end
        end
    end
end

@recipe function f(c,a::Val{:ivc})
    labeling = ["Core","Mantle"]
    coloring = [:blue,:orange]
    for (ic, contact) in enumerate(c.contacts)
        @series begin
            color --> coloring[ic]
            label --> labeling[ic]
            contact
        end
    end
end
@recipe function f(c,a::Val{:coax})
    labeling = ["" for i in 1:19]
    coloring = vcat([:blue], [:orange for i in 1:18])
    for (ic, contact) in enumerate(c.contacts)
        @series begin
            color --> coloring[ic]
            label --> labeling[ic]
            contact
        end
    end
end

@recipe function f(c,a::Val{:segBEGe}, coloring = [])
    isempty(coloring) ? coloring = [:blue, :red, :orange, :grey, :purple] : nothing
    labeling = [ "Core", "Seg. 1", "Seg. 2", "Seg. 3", "Seg. 4"]
    for (ic, contact) in enumerate(c.contacts)
        @series begin
            color --> coloring[ic]
            label --> labeling[ic]
            contact
        end
    end
end

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
