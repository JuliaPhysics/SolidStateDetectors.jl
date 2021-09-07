@recipe function f(t::PointCharge{T}) where {T}

    seriescolor --> :blue
    
    @series begin
        seriescolor --> seriescolor
        seriestype --> :scatter
        markersize --> 6
        label --> "PointCharge"
        t.points
    end
end


@recipe function f(nbc::NBodyChargeCloud{T, N, SH}; connect = true, markersize = 10) where {T, N, SH}
    
    seriescolor --> :blue
    points = nbc.points
    
    #Center point
    @series begin 
        markersize := markersize
        label --> "NBodyChargeCloud"
        seriescolor --> seriescolor
        PointCharge{T}(points[1])
    end
    
    vertex_no = 1
    
    #Shells
    for shell in Base.OneTo(N)
        vertices = get_vertices(nbc.shell_structure)
        @series begin
            markersize := markersize * exp(-(shell))
            label --> ""
            connect --> connect
            seriescolor --> seriescolor
            nbc.shell_structure(points[vertex_no+1:vertex_no+vertices])
        end
        vertex_no += vertices
    end
end

@recipe function f(nbc::NBodyChargeCloud{T, N, NTuple{N, Int}}; connect = true, markersize = 10) where {T, N}
    
    seriescolor --> :blue
    points = nbc.points
    seriestype := :scatter
    
    @series begin 
        markersize := markersize
        label --> "NBodyChargeCloud"
        seriescolor --> seriescolor
        PointCharge{T}(points[1])
    end
    
    vertex_no = 1
    
    #Shells
    for shell in Base.OneTo(N)
        vertices = nbc.shell_structure[shell]
        @series begin
            markersize := markersize * exp(-(shell))
            label --> ""
            connect --> connect
            seriescolor --> seriescolor
            points[vertex_no+1:vertex_no+vertices]
        end
        vertex_no += vertices
    end
end


@recipe function f(t::Tetrahedron{T}; connect = true) where {T}
    
    linestyle --> :dot
    seriescolor --> :blue
    
    @series begin
        seriescolor --> seriescolor
        seriestype --> :scatter
        markersize --> 6
        label --> "Tetrahedron"
        t.points
    end
    
   
    if connect
        @series begin
            seriescolor --> seriescolor
            seriestype --> :line3d
            label --> ""
            linewidth --> 1
            t.points[vcat(1,2,3,4,2)]
        end

        @series begin
            seriescolor --> seriescolor
            seriestype --> :line3d
            label --> ""
            linewidth --> 1
            t.points[vcat(3,1,4)]
        end
    end
end


@recipe function f(i::Octahedron{T}; connect = true) where {T}
    
    linestyle --> :dot
    seriescolor --> :blue
    
    @series begin
        seriescolor --> seriescolor
        seriestype --> :scatter
        markersize --> 6
        label --> "Octahedron"
        i.points
    end
    
    if connect    
        @series begin
            seriescolor --> seriescolor
            seriestype --> :line3d
            label --> ""
            linewidth --> 1
            i.points[vcat(1,2,3,4,5,2,6,5,1,3,6,4,1)]
        end
    end
end


@recipe function f(i::Hexahedron{T}; connect = true) where {T}
    
    linestyle --> :dot
    seriescolor --> :blue
    
    @series begin
        seriescolor --> seriescolor
        seriestype --> :scatter
        markersize --> 6
        label --> "Hexahedron"
        i.points
    end
    
   
    if connect
        @series begin
            seriescolor --> seriescolor
            seriestype --> :line3d
            label --> ""
            linewidth --> 1
            i.points[vcat(1,2,5,3,6,4,7,8,5)]
        end

        @series begin
            seriescolor --> seriescolor
            seriestype --> :line3d
            label --> ""
            linewidth --> 1
            i.points[vcat(3,1,4)]
        end

        @series begin
            seriescolor --> seriescolor
            seriestype --> :line3d
            label --> ""
            linewidth --> 1
            i.points[vcat(2,7)]
        end

        @series begin
            seriescolor --> seriescolor
            seriestype --> :line3d
            label --> ""
            linewidth --> 1
            i.points[vcat(6,8)]
        end
    end
end


@recipe function f(i::Icosahedron{T}; connect = true) where {T}
    
    linestyle --> :dot
    seriescolor --> :blue
    
    @series begin
        seriescolor --> seriescolor
        seriestype --> :scatter
        markersize --> 6
        label --> "Icosahedron"
        i.points
    end
    
   
    if connect
        @series begin
            seriescolor --> seriescolor
            seriestype --> :line3d
            label --> ""
            linewidth --> 1
            i.points[vcat(1,2,7,8,9,10,11,7,3,8,4,9,5,10,6,11,2,3,4,5,6,2)]
        end

        @series begin
            seriescolor --> seriescolor
            seriestype --> :line3d
            label --> ""
            linewidth --> 1
            i.points[vcat(3,1,6)]
        end

        @series begin
            seriescolor --> seriescolor
            seriestype --> :line3d
            label --> ""
            linewidth --> 1
            i.points[vcat(5,1,4)]
        end

        @series begin
            seriescolor --> seriescolor
            seriestype --> :line3d
            label --> ""
            linewidth --> 1
            i.points[vcat(7,12,8)]
        end
        
        @series begin
            seriescolor --> seriescolor
            seriestype --> :line3d
            label --> ""
            linewidth --> 1
            i.points[vcat(9,12,10)]
        end
        
        @series begin
            seriescolor --> seriescolor
            seriestype --> :line3d
            label --> ""
            linewidth --> 1
            i.points[vcat(11,12)]
        end
    end
end


@recipe function f(i::Dodecahedron{T}; connect = true) where {T}
    
    linestyle --> :dot
    seriescolor --> :blue
    
    @series begin
        seriescolor --> seriescolor
        seriestype --> :scatter
        markersize --> 6
        label --> "Dodecahedron"
        i.points
    end
    
   
    if connect
        @series begin
            seriescolor --> seriescolor
            seriestype --> :line3d
            label --> ""
            linewidth --> 1
            i.points[vcat(1,2,3,4,5,1,6,11,7,12,8,13,9,14,10,15,6)]
        end

        @series begin
            seriescolor --> seriescolor
            seriestype --> :line3d
            label --> ""
            linewidth --> 1
            i.points[vcat(15,20,16,17,18,19,20)]
        end
        
        for j in [2,3,4,5,11,12,13,14]
            @series begin
                seriescolor --> seriescolor
                seriestype --> :line3d
                label --> ""
                linewidth --> 1
                i.points[vcat(j,j+5)]
            end
        end
    end
end