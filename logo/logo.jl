#Versions:
#Plots v1.40.2
#PyPlot v2.11.2
using Plots; pyplot()

#colors extracted from the Julia logo
jred = [RGB(213/255,99/255,92/255), RGB(203/255,60/255,51/255)]
jgreen = [RGB(96/255,173/255,81/255), RGB(56/255,152/255,38/255)]
jblue = [RGB(102/255,130/255, 233/255), RGB(64/255,99/255,216/255)]
jpurple = [RGB(170/255,121/255,193/255),RGB(149/255,88/255,178/255)]
jtext = [[jpurple[2],jpurple[2],jpurple[2]],[jpurple[2],jpurple[2],jpurple[2]]]

function plot_frame(bc = :transparent, negative_color = :white)

    #detector geometry
    top=3.5
    top2=2
    bot=6.8
    inner=2
    middle=5.3
    outer=6.2
    ratio=4
    lwidth=5

    plot(xlims = (-6.5,8.5), ylims = (-9.15,5.85), legend = false, size = (429,429),
    showaxis = false, grid = false, background_color = bc)
    if bc == :transparent
        #plot detector background
        rectangle(w, h, x, y) = Shape(x .+ [0,w,w,0], y .+ [0,0,h,h])
        plot!(rectangle(2*outer,top2+bot,-outer,-bot), color = negative_color, linecolor = :transparent)
        plot!(Shape([-outer,outer,middle,-middle],[top2-0.1,top2-0.1,top,top]), color = negative_color, linecolor = :transparent)
        for i in inner:0.1:middle; plot!([i*cos(t) for t in 0:pi/100:3pi], [i/ratio*sin(t) + top for t in 0:pi/100:3pi], linewidth = 3, color = negative_color); end
        for i in 0:0.1:outer; plot!([i*cos(t) for t in 0:pi/100:3pi], [i/ratio*sin(t) - bot for t in 0:pi/100:3pi], linewidth = 3, color = negative_color); end
    end
    plot!()
end

function plot_detector(negative_color = :white)

    #detector geometry
    top=3.5
    top2=2
    bot=6.8
    inner=2
    middle=5.3
    outer=6.2
    ratio=4
    lwidth=5

    #plot detector lines and contact
    for i in 0:0.05:0.95 plot!([i*inner*cos(t) for t in 0:pi/100:3pi], [i*inner/ratio*sin(t) + top for t in 0:pi/100:3pi], linewidth = 3, color = jblue[1]) end
    plot!([inner*cos(t) for t in 0:pi/100:3pi], [inner/ratio*sin(t)+top for t in 0:pi/100:3pi], linewidth = lwidth, color = jblue[2])
    plot!([middle*cos(t) for t in 0:pi/100:3pi], [middle/ratio*sin(t)+top for t in 0:pi/100:3pi], linewidth = lwidth, color = :gray)
    plot!(vcat([-middle,-outer],[outer*cos(t) for t in -pi:pi/100:0],[outer,middle]), vcat([top,top2],[outer/ratio*sin(t)-bot for t in -pi:pi/100:0],[top2,top]), linewidth = lwidth, color = :gray)
end


function SSD(k::Integer)
    if k >= 420 return 3 end
    if k >= 300 return 2 end
    if k >= 150 return 1 end
end

function plot_rest(k::Integer, label::Bool, animate::Bool, negative_color = :white)

    #detector geometry
    top=3.5
    top2=2
    bot=6.8
    inner=2
    middle=5.3
    outer=6.2
    ratio=4

    #define incoming γ ray
    γmax = 9      #length of γ ray
    ϕ0 = -0.8-pi  #phase of γ ray
    α= -pi*1/2.8  #rotation angle for γ ray
    x = collect(γmax:-0.02:2.9)
    y = sin.(x * 0.7pi .+ ϕ0)
    newx = vcat(4.3 .+ x * cos(α) .+ y * sin(α), [4.8+1.9*cos(t) for t in pi/2-0.1:-pi/100:-pi/1.5])
    newy = vcat(-4.55 .- x * sin(α) .+ y * cos(α), [-3.25+1.6*sin(t) for t in pi/2-0.1:-pi/100:-pi/1.5])

    lwidth = 5    #width of lines

    #plot drift paths of electrons and holes
    if k > 400
        kk = k - 400
        hx = [outer*cos(t) for t in -2pi/3.2:-pi/300:-pi][1:min(kk,end)]
        hy = [outer/4*sin(t)-3.6 for t in -2pi/3.2:-pi/300:-pi][1:min(kk,end)]
        ex = [0,0]
        ey = [sqrt(25-6.25)-5, min(sqrt(25-6.25)-5 + (top - sqrt(25-6.25)+5)*kk/length(-2pi/3.2:-pi/300:-pi),top)]
        plot!(ex,ey, width = lwidth, color = jgreen[2])
        plot!(hx,hy, width = lwidth, color = jred[2])
    end
    plot_detector(negative_color)

    #plot incoming γ ray
    contourx = newx[1:min(k,end)]
    contoury = newy[1:min(k,end)]
    idx = findall(.!(outer - 0.127 .< contourx .< outer + 0.13))
    for i in findall(diff(idx) .> 1)
        line_idx = idx[i]:idx[i+1]
        dymax = abs(0.2*sqrt(1 + ((-)(extrema(contoury[line_idx])...)/(-)(extrema(contourx[line_idx])...))^2))
        for dy in range(-dymax, stop = dymax, length = ceil(Int,dymax/0.001))
            plot!(contourx[line_idx], contoury[line_idx] .+ dy, width = 0.2*lwidth, color = negative_color)
        end
    end
    if animate
        scatter!(newx[1:min(k,end)], newy[1:min(k,end)], markersize = 10, markerstrokewidth = 0, color = negative_color)
    end
    plot!(newx[1:min(k,end)], newy[1:min(k,end)], width = lwidth, color = jtext[2][2])
    


    #plot the three Julia dots
    if k > 350
        msize = min(k - 350,78)
        scatter!([2.5],[-5], color = jpurple[1], markerstrokecolor = jpurple[2], markerstrokewidth = lwidth, markersize = msize)
        scatter!([-2.5],[-5], color = jred[1], markerstrokecolor = jred[2], markerstrokewidth = lwidth, markersize = msize)
        scatter!([0],[sqrt(25-6.25)-5], color = jgreen[1], markerstrokecolor = jgreen[2], markerstrokewidth = lwidth, markersize = msize)
    end

    #add labels
    if label
        if k >= 150
            annotate!([(8.3, 1.57, Plots.text("olid", 50, jtext[1][1], :left)),
               (7.05,-1.05, Plots.text("tate", 50, jtext[1][2], :left)),
               (6.8,-3.7, Plots.text("etectors", 50, jtext[1][3], :left))][1:SSD(k)])
        end

        plot!(xlims = (-6.5,16.7), ylims = (-8.5,5.2), size = (662,391))
    end
    plot!()
end


function get_logo(; animate::Bool=false, label::Bool=true, transparent::Bool = true, dark::Bool = false, step::Int = 5)
    negative_color = dark ? :black : :white
    if animate
        ks::Vector{Int} = collect(0:step:600)
        if !(600 in ks) push!(ks, 600) end
        return @animate for k in ks
            k = (k == 0 ? k = 600 : k)
            plot_frame(negative_color)
            plot_rest(k,label,animate,negative_color)
        end
    else
        plot_frame(transparent ? :transparent : negative_color, negative_color)
        plot_rest(600,label,animate,negative_color)
    end
end


"""
    get_and_save_logo( [; animate::Bool=false, label::Bool=true])

Compute and save the SolidStateDetectors.jl logo.

# Keyword arguments
- `animate::Bool=false`: if set to `true`, the logo will be saved as animated GIF, otherwise as static SVG file.
- `label::Bool=true`: if set to `true`, the words "Solid State Detectors" will appear on the right of the logo.
- `dark::Bool = false`: if set to `true`, the logo cutouts and background will be plotted black (for dark mode).
- `transparent::Bool=true`: if set to `true`, the background will be plotted transparent

"""
function get_and_save_logo(; animate::Bool = false, label::Bool = true, transparent::Bool = true, dark::Bool = false, fps = 15, step::Int = 5, loop::Int = 0)
    logo = get_logo(animate = animate, label = label, transparent = transparent, dark = dark, step = step);
    basename = "logo" * (dark ? "_dark" : "") * (label ? "" : "_no_name")
    animate ? gif(logo, basename*".gif", fps = fps, loop = loop) : savefig(basename*".svg")
end

# get_and_save_logo(animate = true, label = false, fps = 300, step = 5, loop = false)
