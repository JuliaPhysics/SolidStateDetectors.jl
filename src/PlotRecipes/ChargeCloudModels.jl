@recipe function f(t::PointCharge{T}) where {T}

    seriescolor --> :blue
    
    @series begin
        seriescolor --> seriescolor
        seriestype --> :scatter
        markersize --> 6
        label --> "PointCharge"
        t.locations
    end
end


@recipe function f(nbc::NBodyChargeCloud{T, N, SH}; connect = true, markersize = 10) where {T, N, SH}
    
    seriescolor --> :blue
    points = nbc.locations
    
    #Center point
    @series begin 
        markersize := markersize
        label --> "NBodyChargeCloud"
        seriescolor --> seriescolor
        PointCharge(points[1])
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
    points = nbc.locations
    seriestype := :scatter
    
    @series begin 
        markersize := markersize
        label --> "NBodyChargeCloud"
        seriescolor --> seriescolor
        PointCharge(points[1])
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
        t.locations
    end
    
   
    if connect
        @series begin
            seriescolor --> seriescolor
            seriestype --> :line3d
            label --> ""
            linewidth --> 1
            t.locations[vcat(1,2,3,4,2)]
        end

        @series begin
            seriescolor --> seriescolor
            seriestype --> :line3d
            label --> ""
            linewidth --> 1
            t.locations[vcat(3,1,4)]
        end
    end
end


@recipe function f(o::Octahedron{T}; connect = true) where {T}
    
    linestyle --> :dot
    seriescolor --> :blue
    
    @series begin
        seriescolor --> seriescolor
        seriestype --> :scatter
        markersize --> 6
        label --> "Octahedron"
        o.locations
    end
    
    if connect    
        @series begin
            seriescolor --> seriescolor
            seriestype --> :line3d
            label --> ""
            linewidth --> 1
            o.locations[vcat(1,2,3,4,5,2,6,5,1,3,6,4,1)]
        end
    end
end


@recipe function f(h::Hexahedron{T}; connect = true) where {T}
    
    linestyle --> :dot
    seriescolor --> :blue
    
    @series begin
        seriescolor --> seriescolor
        seriestype --> :scatter
        markersize --> 6
        label --> "Hexahedron"
        h.locations
    end
    
   
    if connect
        @series begin
            seriescolor --> seriescolor
            seriestype --> :line3d
            label --> ""
            linewidth --> 1
            h.locations[vcat(1,2,5,3,6,4,7,8,5)]
        end

        @series begin
            seriescolor --> seriescolor
            seriestype --> :line3d
            label --> ""
            linewidth --> 1
            h.locations[vcat(3,1,4)]
        end

        @series begin
            seriescolor --> seriescolor
            seriestype --> :line3d
            label --> ""
            linewidth --> 1
            h.locations[vcat(2,7)]
        end

        @series begin
            seriescolor --> seriescolor
            seriestype --> :line3d
            label --> ""
            linewidth --> 1
            h.locations[vcat(6,8)]
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
        i.locations
    end
    
   
    if connect
        @series begin
            seriescolor --> seriescolor
            seriestype --> :line3d
            label --> ""
            linewidth --> 1
            i.locations[vcat(1,2,7,8,9,10,11,7,3,8,4,9,5,10,6,11,2,3,4,5,6,2)]
        end

        @series begin
            seriescolor --> seriescolor
            seriestype --> :line3d
            label --> ""
            linewidth --> 1
            i.locations[vcat(3,1,6)]
        end

        @series begin
            seriescolor --> seriescolor
            seriestype --> :line3d
            label --> ""
            linewidth --> 1
            i.locations[vcat(5,1,4)]
        end

        @series begin
            seriescolor --> seriescolor
            seriestype --> :line3d
            label --> ""
            linewidth --> 1
            i.locations[vcat(7,12,8)]
        end
        
        @series begin
            seriescolor --> seriescolor
            seriestype --> :line3d
            label --> ""
            linewidth --> 1
            i.locations[vcat(9,12,10)]
        end
        
        @series begin
            seriescolor --> seriescolor
            seriestype --> :line3d
            label --> ""
            linewidth --> 1
            i.locations[vcat(11,12)]
        end
    end
end


@recipe function f(d::Dodecahedron{T}; connect = true) where {T}
    
    linestyle --> :dot
    seriescolor --> :blue
    
    @series begin
        seriescolor --> seriescolor
        seriestype --> :scatter
        markersize --> 6
        label --> "Dodecahedron"
        d.locations
    end
    
   
    if connect
        @series begin
            seriescolor --> seriescolor
            seriestype --> :line3d
            label --> ""
            linewidth --> 1
            d.locations[vcat(1,2,3,4,5,1,6,11,7,12,8,13,9,14,10,15,6)]
        end

        @series begin
            seriescolor --> seriescolor
            seriestype --> :line3d
            label --> ""
            linewidth --> 1
            d.locations[vcat(15,20,16,17,18,19,20)]
        end
        
        for j in [2,3,4,5,11,12,13,14]
            @series begin
                seriescolor --> seriescolor
                seriestype --> :line3d
                label --> ""
                linewidth --> 1
                d.locations[vcat(j,j+5)]
            end
        end
    end
end

@recipe function f(m::AbstractParticleSource; length = 0.01)
    
    if hasproperty(m, :opening_angle)
        if iszero(m.opening_angle)
            v = m.position
            vnew = m.position + length * normalize(m.direction)
            @series begin
                linewidth --> 2
                color --> :green
                label := ""
                [v.x, vnew.x], [v.y, vnew.y], [v.z, vnew.z]
            end
        elseif m.opening_angle <= 90u"°"
            d = normalize(m.direction)
            a = normalize(d × (abs(d.x) == 1 ? CartesianVector(0,1,0) : CartesianVector(1,0,0)))
            b = normalize(a × d)
            rot = hcat(a,b,d)
            cone = SolidStateDetectors.ConstructiveSolidGeometry.Cone(r = ((0,0),(0,length*sin(m.opening_angle))), hZ = length*cos(m.opening_angle)/2, 
            origin = rot * [0,0,length*cos(m.opening_angle)/2] + m.position, 
            rotation = rot)

            @series begin
                linewidth := 0
                color --> :green
                label := ""
                cone
            end
        end
    end
    
    @series begin
        seriescolor := :gray
        markersize --> 5
        markerstrokewidth --> 1
        linecolor := :black
        label --> string(typeof(m).name.name)
        m.position
    end
end