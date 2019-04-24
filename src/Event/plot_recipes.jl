
# @recipe function f(e::Event;coloring=[], labeling=[],show_pulses=true)
#     if show_pulses==true
#         n_wps = size(e.weighting_potentials,1)
#         width = 800
#         length = width+n_wps*width/2
#         size  -->  (width,length)
#         myheights = [width/length]
#         for i in 1:n_wps
#             push!(myheights,(width/2)/length)
#         end

#         layout  -->  (n_wps+1,1) #grid(n_wps+1,1,heights=myheights)
#         @series begin
#             subplot := 1
#             coloring --> coloring
#             labeling --> labeling
#             e.detector
#         end
#         for itr in 1:e.n_sites
#             @series begin
#                 itr == 1 ? showlabel --> true : showlabel --> false
#                 subplot := 1
#                 myscale := 1/e.detector.geometry_unit_factor
#                 e.trajectories_e[itr]
#             end
#             @series begin
#                 itr == 1 ? showlabel --> true : showlabel --> false
#                 subplot := 1
#                 myscale := 1/e.detector.geometry_unit_factor
#                 e.trajectories_h[itr]
#             end
#         end
#         for iwp in 1:n_wps
#             @series begin
#                 subplot := 1+iwp
#                 e.pulses[iwp]
#             end
#         end
#     else
#         width,length = 800, 800
#         size  -->  (width,length)
#         @series begin
#             subplot := 1
#             coloring --> coloring
#             labeling --> labeling
#             e.detector
#         end

#         for itr in 1:e.n_sites
#             @series begin
#                 itr == 1 ? showlabel --> true : showlabel --> false
#                 subplot := 1
#                 myscale := 1/e.detector.geometry_unit_factor
#                 e.trajectories_e[itr]
#             end
#             @series begin
#                 itr == 1 ? showlabel --> true : showlabel --> false
#                 subplot := 1
#                 myscale := 1/e.detector.geometry_unit_factor
#                 e.trajectories_h[itr]
#             end
#         end
#     end
# end
