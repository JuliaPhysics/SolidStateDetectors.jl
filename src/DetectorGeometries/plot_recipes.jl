
@recipe function f( d::InvertedCoax, dim::Symbol;
                    r=missing,
                    θ=missing,
                    z=missing,
                    half=true)
    T = typeof(d.crystal_radius)

    cross_section::Symbol, value::T = if ismissing(θ) && ismissing(r) && ismissing(z)
        :θ, 0
    elseif !ismissing(θ) && ismissing(r) && ismissing(z)
        :θ, θ
    elseif ismissing(θ) && !ismissing(r) && ismissing(z)
        :r, r
    elseif ismissing(θ) && ismissing(r) && !ismissing(z)
        :z, z
    else
        error(ArgumentError, ": Only one of the keywords `r, θ, z` is allowed.")
    end

    if cross_section == :θ
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

@recipe function f( d::BEGe, dim::Symbol;
                    r=missing,
                    θ=missing,
                    z=missing,
                    half=true)
    T = typeof(d.crystal_radius)

    cross_section::Symbol, value::T = if ismissing(θ) && ismissing(r) && ismissing(z)
        :θ, 0
    elseif !ismissing(θ) && ismissing(r) && ismissing(z)
        :θ, θ
    elseif ismissing(θ) && !ismissing(r) && ismissing(z)
        :r, r
    elseif ismissing(θ) && ismissing(r) && !ismissing(z)
        :z, z
    else
        error(ArgumentError, ": Only one of the keywords `r, θ, z` is allowed.")
    end

    if cross_section == :θ
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
@recipe function f(c::Coax)
    size:=(800,800)
    linewidth:=2
    linecolor:=:black
    label-->""
    aspect_ratio --> 1
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

struct outer_taper{T<:AbstractFloat}
    rStart::T
    rStop::T
    θStart::T
    θStop::T
    zStart::T
    zStop::T
    orientation::String
    function outer_taper(rStart::T, rStop::T, θStart::T, θStop::T, zStart::T, zStop::T, orientation) where T<:AbstractFloat
        return new{T}(rStart, rStop, θStart, θStop, zStart, zStop, orientation)
    end
end

function outer_taper(rStart, rStop, θStart, θStop, zStart, zStop, orientation)
    return outer_taper{typeof(rStart)}(rStart, rStop, θStart, θStop, zStart, zStop, orientation)
end

@recipe function f(b::InvertedCoax; coloring=[], labeling=[])
    if b.name =="ExampleInvertedCoax"
        coloring = [:blue, :orange, :orange, :orange, :orange, :orange, :orange, :orange]
        labeling = ["Core", "Mantle", "", "", "", "", "", ""]
    end
    aspect_ratio --> 1
    f=1/b.geometry_unit_factor
    for i_part in 1:size(b.segmentation_r_ranges,1)
        rStart = f*b.segmentation_r_ranges[i_part][1]
        rStop = f*b.segmentation_r_ranges[i_part][2]
        θStart = b.segmentation_phi_ranges[i_part][1]
        θStop = b.segmentation_phi_ranges[i_part][2]
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
                Tubs("",1,1.,"",0.0,rStart,rStop,θStart,θStop,zStart,zStop)
            end
        else
            @series begin
                if !isempty(labeling)
                    label := labeling[i_part]
                else
                    label := ""
                end
                # Taper("",1,1.,"",0.0,rStart,rStop,θStart,θStop,zStart,zStop,"bl")
                outer_taper(rStart,rStop,θStart,θStop,zStart,zStop,"c//")
            end
        end
    end
end

@recipe function f(b::BEGe; coloring=[], labeling=[])
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
    aspect_ratio --> 1
    f=1/b.geometry_unit_factor
    for i_part in 1:size(b.segmentation_r_ranges,1)
        rStart = f*b.segmentation_r_ranges[i_part][1]
        rStop = f*b.segmentation_r_ranges[i_part][2]
        θStart = b.segmentation_phi_ranges[i_part][1]
        θStop = b.segmentation_phi_ranges[i_part][2]
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
                Tubs("",1,1.,"",0.0,rStart,rStop,θStart,θStop,zStart,zStop)
            end
        else
            @series begin
                if !isempty(labeling)
                    label := labeling[i_part]
                else
                    label := ""
                end
                # Taper("",1,1.,"",0.0,rStart,rStop,θStart,θStop,zStart,zStop,"bl")
                outer_taper(rStart,rStop,θStart,θStop,zStart,zStop,"c//")
            end
        end
    end
end

@recipe function f(Vol::Tubs{T}) where T<:AbstractFloat

    @series begin
        partialcircle_3d(Vol.rStop,Vol.θStart,Vol.θStop,[0,0,Vol.zStart])
    end
    @series begin
        label:= ""
        partialcircle_3d(Vol.rStop,Vol.θStart,Vol.θStop,[0,0,Vol.zStop])
    end
    @series begin
        label:= ""
        partialcircle_3d(Vol.rStart,Vol.θStart,Vol.θStop,[0,0,Vol.zStart])
    end
    @series begin
        label:= ""
        partialcircle_3d(Vol.rStart,Vol.θStart,Vol.θStop,[0,0,Vol.zStop])
    end
    ## Vertical Lines

    if !iszero(Vol.rStart)
        @series begin
            label:= ""
            line_3d(Vol.rStart,Vol.rStart,Vol.θStart,Vol.θStart,Vol.zStart,Vol.zStop)
        end
    end
    if !iszero(Vol.rStart)
        @series begin
            label:= ""
            line_3d(Vol.rStart,Vol.rStart,Vol.θStop,Vol.θStop,Vol.zStart,Vol.zStop)
        end
    end
    @series begin
        label:= ""
        line_3d(Vol.rStop,Vol.rStop,Vol.θStart,Vol.θStart,Vol.zStart,Vol.zStop)
    end
    @series begin
        label:= ""
        line_3d(Vol.rStop,Vol.rStop,Vol.θStop,Vol.θStop,Vol.zStart,Vol.zStop)
    end

    ##Horizontal Lines

    if !isapprox((Vol.θStop - Vol.θStart)%2π , 0.0,atol=0.00001)
        @series begin
            label:= ""
            line_3d(Vol.rStart,Vol.rStop,Vol.θStart,Vol.θStart,Vol.zStart,Vol.zStart)
        end
        @series begin
            label:= ""
            line_3d(Vol.rStart,Vol.rStop,Vol.θStop,Vol.θStop,Vol.zStart,Vol.zStart)
        end
        @series begin
            label:= ""
            line_3d(Vol.rStart,Vol.rStop,Vol.θStart,Vol.θStart,Vol.zStop,Vol.zStop)
        end
        @series begin
            label:= ""
            line_3d(Vol.rStart,Vol.rStop,Vol.θStop,Vol.θStop,Vol.zStop,Vol.zStop)
        end
    end
end

@recipe function f(o::outer_taper,n_aux_lines =0)
    if o.orientation=="c//"
        @series begin
            line_3d(o.rStop,o.rStart,o.θStart,o.θStart,o.zStart,o.zStop)
        end
        @series begin
            line_3d(o.rStop,o.rStart,o.θStop,o.θStop,o.zStart,o.zStop)
        end
        for ia in 0:n_aux_lines
            @series begin
                line_3d(o.rStop,o.rStart,o.θStart+ia*(o.θStop-o.θStart)/(n_aux_lines+1),o.θStart+ia*(o.θStop-o.θStart)/(n_aux_lines+1),o.zStart,o.zStop)
            end
        end
    end
end

function connect_points()
    nothing
end

function mylinspace(Start,Stop,nSteps)
    return collect(Start:(Stop-Start)/nSteps:Stop)
end

function partialcircle_3d(radius,phiStart,phiStop,Translate::Vector;nSteps=400)
    phirange = mylinspace(phiStart,phiStop,nSteps)

    x::Vector{AbstractFloat}=map(x->radius*cos.(x),phirange)
    y::Vector{AbstractFloat}=map(x->radius*sin.(x),phirange)
    # z::Vector{AbstractFloat}=[Translate[3] for i in 1:nSteps]
    z::Vector{AbstractFloat}=map(x->Translate[3],phirange)
    x.+=Translate[1]
    y.+=Translate[2]

    return x,y,z
end

function partialcircle(radius,phiStart,phiStop,Translate::Vector;nSteps=400)
    phirange = mylinspace(phiStart,phiStop,nSteps)
    x::Vector{AbstractFloat}=map(x->radius*cos.(x),phirange)
    y::Vector{AbstractFloat}=map(x->radius*sin.(x),phirange)
    x.+=Translate[1]
    y.+=Translate[2]
    return x,y
end

function line_3d(r1,r2,phi1,phi2,z1,z2)
    x1::Vector = [r1*cos(phi1),r1*sin(phi1),z1]
    x2::Vector = [r2*cos(phi2),r2*sin(phi2),z2]
    return [x1[1],x2[1]] , [x1[2],x2[2]], [x1[3],x2[3]]
end
