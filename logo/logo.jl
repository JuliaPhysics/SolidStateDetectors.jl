#Versions:
#Plots v0.28.2
#PyPlot v2.8.2
using Plots; pyplot()

#colors extracted from the Julia logo
jred = [RGB(213/255,99/255,92/255), RGB(203/255,60/255,51/255)]
jgreen = [RGB(96/255,173/255,81/255), RGB(56/255,152/255,38/255)]
jblue = [RGB(102/255,130/255, 233/255), RGB(64/255,99/255,216/255)]
jpurple = [RGB(170/255,121/255,193/255),RGB(149/255,88/255,178/255)]
jtext = [[jpurple[2],jpurple[2],jpurple[2]],[jpurple[2],jpurple[2],jpurple[2]]]
jdet = [:white, :gray]

function plot_frame(bc=:transparent)

    #detector geometry
    top=3.5
    top2=2
    bot=6.8
    inner=2
    middle=5.3
    outer=6.3
    ratio=4
    lwidth=5

    plot(xlims = (-6.5,8.5), ylims = (-9.15,5.85), legend = false, size = (429,429),
    showaxis = false, grid = false, background_color = bc)
    if bc == :transparent
        #plot detector background
        rectangle(w, h, x, y) = Shape(x .+ [0,w,w,0], y .+ [0,0,h,h])
        plot!(rectangle(2*outer,top2+bot,-outer,-bot), color = jdet[1], linecolor = :transparent)
        plot!(Shape([-outer,outer,middle,-middle],[top2-0.1,top2-0.1,top,top]), color = jdet[1], linecolor = :transparent)
        for i in inner:0.1:middle; plot!([i*cos(t) for t in 0:pi/100:2pi], [i/ratio*sin(t) + top for t in 0:pi/100:2pi], linewidth = 3, color = jdet[1]); end
        for i in 0:0.1:outer; plot!([i*cos(t) for t in 0:pi/100:2pi], [i/ratio*sin(t) - bot for t in 0:pi/100:2pi], linewidth = 3, color = jdet[1]); end
    end
    plot!()
end

function plot_detector()

    #detector geometry
    top=3.5
    top2=2
    bot=6.8
    inner=2
    middle=5.3
    outer=6.3
    ratio=4
    lwidth=5

    #plot detector lines and contact
    for i in 0:0.05:0.95 plot!([i*inner*cos(t) for t in 0:pi/100:2pi], [i*inner/ratio*sin(t) + top for t in 0:pi/100:2pi], linewidth = 3, color = jblue[1]) end
    plot!([inner*cos(t) for t in 0:pi/100:2pi], [inner/ratio*sin(t)+top for t in 0:pi/100:2pi], linewidth = lwidth, color = jblue[2])
    plot!([middle*cos(t) for t in 0:pi/100:2pi], [middle/ratio*sin(t)+top for t in 0:pi/100:2pi], linewidth = lwidth, color = jdet[2])
    plot!(vcat([-middle,-outer],[outer*cos(t) for t in -pi:pi/100:0],[outer,middle]), vcat([top,top2],[outer/ratio*sin(t)-bot for t in -pi:pi/100:0],[top2,top]), linewidth = lwidth, color = jdet[2])
end


function SSD(k::Integer)
    if k >= 420 return 3 end
    if k >= 300 return 2 end
    if k >= 150 return 1 end
end

function plot_rest(k::Integer, label::Bool)

    #detector geometry
    top=3.5
    top2=2
    bot=6.8
    inner=2
    middle=5.3
    outer=6.3
    ratio=4

    #plot drift paths of electrons and holes
    ex = [0,0]
    ey = [sqrt(25-6.25)-5, top]
    hx = [6.3*cos(t) for t in -2pi/3.2:-pi/100:-pi]
    hy = [6.3/4*sin(t)-3.6 for t in -2pi/3.2:-pi/100:-pi]

    #plot incoming γ ray
    γmax = 9      #length of γ ray
    ϕ0 = -0.8-pi  #phase of γ ray
    α= -pi*1/2.8  #rotation angle for γ ray
    x = collect(γmax:-0.02:2.9)
    y = sin.(x * 0.7pi .+ ϕ0)
    newx = vcat(4 .+ x * cos(α) .+ y * sin(α), [4.5+1.5*cos(t) for t in pi/2-0.1:-pi/100:-pi/1.5])
    newy = vcat(-4 .- x * sin(α) .+ y * cos(α), [-2.6+1.5*sin(t) for t in pi/2-0.1:-pi/100:-pi/1.5])

    lwidth = 5    #width of lines

    if k > 400
        kk = k - 400
        hx = [6.3*cos(t) for t in -2pi/3.2:-pi/300:-pi][1:min(kk,end)]
        hy = [6.3/4*sin(t)-3.6 for t in -2pi/3.2:-pi/300:-pi][1:min(kk,end)]
        ex = [0,0]
        ey = [sqrt(25-6.25)-5, min(sqrt(25-6.25)-5 + (top - sqrt(25-6.25)+5)*kk/length(-2pi/3.2:-pi/300:-pi),top)]
        plot!(ex,ey, width = lwidth, color = jgreen[2])
        plot!(hx,hy, width = lwidth, color = jred[2])
    end
    plot_detector()

    plot!(newx[1:min(k,length(x))], newy[1:min(k,length(x))], width = 2*lwidth, color = :white)
    plot!(newx[1:min(k,end)], newy[1:min(k,end)], width = lwidth, color = jtext[2][2])


    #plot the three Julia dots
    if k > 350
        msize = min(k - 350,78)
        scatter!([2.5],[-5], color = jpurple[1], markerstrokecolor = jpurple[2], markerstrokewidth = lwidth, markersize = msize)
        scatter!([-2.5],[-5], color = jred[1], markerstrokecolor = jred[2], markerstrokewidth = lwidth, markersize = msize)
        scatter!([0],[sqrt(25-6.25)-5], color = jgreen[1], markerstrokecolor = jgreen[2], markerstrokewidth = lwidth, markersize = msize)
    end

    if label
        if k >= 150
            annotate!([(8.15, 2.1, Plots.text("olid", 50, jtext[1][1], :left)),
               (6.9,-0.45, Plots.text("tate", 50, jtext[1][2], :left)),
               (6.4,-3.1, Plots.text("etectors", 50, jtext[1][3], :left))][1:SSD(k)])
        end

        plot!(xlims = (-6.5,16.2), ylims = (-8.5,5.2), size = (649,391))
    end
    plot!()
end

function get_logo(; animate::Bool=false, label::Bool=true)
    if animate
        return @animate for k in 1:600
            plot_frame(:white)
            plot_rest(k, label)
        end
    else
        plot_frame()
        plot_rest(600,label)
    end
end


"""
    get_and_save_logo( [; animate::Bool=false, label::Bool=true])

Compute and save the SolidStateDetectors.jl logo.

# Keyword arguments
- `animate::Bool=false`: if set to `true`, the logo will be saved as animated GIF, otherwise as static SVG file.
- `label::Bool=true`: if set to `true`, the words "Solid State Detectors" will appear on the right of the logo.

"""
function get_and_save_logo(; animate::Bool=false, label::Bool=true)
    logo = get_logo(animate = animate, label = label);
    basename = "logo" * (label ? "" : "_no_name")
    animate ? gif(logo, basename*".gif") : savefig(basename*".svg")
end

get_and_save_logo()
