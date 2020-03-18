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

#jtext = [[:black,:black,:black],[:black,:black,:black]]
#jtext = [[jgreen[1],jred[1],jpurple[1]],[jgreen[2],jred[2],jpurple[2]]]
#jdet = [:lightgray, :gray]

plot(xlims = (-6.5,8.5), ylims = (-9.15,5.85), legend = false, size = (429,429),
    showaxis = false, grid = false, background_color = :transparent)

msize = 78    #size of dots
lwidth = 5    #width of lines
γmax = 9      #length of γ ray
ϕ0 = -0.8-pi  #phase of γ ray
α= -pi*1/2.8  #rotation angle for γ ray
label = true  #show the text "Solid State Detectors"

#detector geometry
top=3.5
top2=2
bot=6.8
inner=2
middle=5.3
outer=6.3
ratio=4


#plot detector background
rectangle(w, h, x, y) = Shape(x .+ [0,w,w,0], y .+ [0,0,h,h])
plot!(rectangle(2*outer,top2+bot,-outer,-bot), color = jdet[1], linecolor = :transparent)
plot!(Shape([-outer,outer,middle,-middle],[top2-0.1,top2-0.1,top,top]), color = jdet[1], linecolor = :transparent)
for i in inner:0.1:middle; plot!([i*cos(t) for t in 0:pi/100:2pi], [i/ratio*sin(t) + top for t in 0:pi/100:2pi], linewidth = 3, color = jdet[1]); end
for i in 0:0.1:outer; plot!([i*cos(t) for t in 0:pi/100:2pi], [i/ratio*sin(t) - bot for t in 0:pi/100:2pi], linewidth = 3, color = jdet[1]); end


#plot drift paths of electrons and holes
plot!([0,0],[sqrt(25-6.25)-5, top], width = lwidth, color = jgreen[2])
plot!([6.3*cos(t) for t in -pi:pi/100:-2pi/3],[6.3/4*sin(t)-3.6 for t in -pi:pi/100:-2pi/3], width = lwidth, color = jred[2])
plot!([4.5+1.5*cos(t) for t in pi/2+0.1:-pi/100:-pi/1.5],[-2.6+1.5*sin(t) for t in pi/2+0.1:-pi/100:-pi/1.5], width=lwidth, color = jtext[2][3])

#plot the three Julia dots
scatter!([2.5],[-5], color = jpurple[1], markerstrokecolor = jpurple[2], markerstrokewidth = lwidth, markersize = msize)
scatter!([-2.5],[-5], color = jred[1], markerstrokecolor = jred[2], markerstrokewidth = lwidth, markersize = msize)
scatter!([0],[sqrt(25-6.25)-5], color = jgreen[1], markerstrokecolor = jgreen[2], markerstrokewidth = lwidth, markersize = msize)

#plot detector lines and contact
for i in 0:0.05:0.95 plot!([i*inner*cos(t) for t in 0:pi/100:2pi], [i*inner/ratio*sin(t) + top for t in 0:pi/100:2pi], linewidth = 3, color = jblue[1]) end
plot!([inner*cos(t) for t in 0:pi/100:2pi], [inner/ratio*sin(t)+top for t in 0:pi/100:2pi], linewidth = lwidth, color = jblue[2])
plot!([middle*cos(t) for t in 0:pi/100:2pi], [middle/ratio*sin(t)+top for t in 0:pi/100:2pi], linewidth = lwidth, color = jdet[2])
plot!(vcat([-middle,-outer],[outer*cos(t) for t in -pi:pi/100:0],[outer,middle]), vcat([top,top2],[outer/ratio*sin(t)-bot for t in -pi:pi/100:0],[top2,top]), linewidth = lwidth, color = jdet[2])

#plot the incoming γ ray
x = collect(2.95:0.02:γmax)
y = sin.(x * 0.7pi .+ ϕ0)
newx = 4 .+ x * cos(α) .+ y * sin(α)
newy = -4 .- x * sin(α) .+ y * cos(α)
plot!(newx, newy, width = 2*lwidth, color = :white)
plot!(newx, newy, width = lwidth, color = jtext[2][2])
plot!(newx[round(Int,end*0.53):end], newy[round(Int,end*0.53):end], width = lwidth, color = jtext[2][1])
plot!([4.5+1.5*cos(t) for t in pi/2+0.1:-pi/100:0],[-2.6+1.5*sin(t) for t in pi/2+0.1:-pi/100:0], width=lwidth, color = jtext[2][3])

#add text and adjust size of the logo
if label
annotate!([(8.15, 2.1, Plots.text("olid", 50, jtext[1][1], :left)),
           (6.9,-0.45, Plots.text("tate", 50, jtext[1][2], :left)),
           (6.3,-3.1, Plots.text("etectors", 50, jtext[1][3], :left))])
plot!(xlims = (-6.5,16.2), ylims = (-8.5,5.2), size = (649,391), showaxis = false)
end

savefig("ssd_logo.svg")
