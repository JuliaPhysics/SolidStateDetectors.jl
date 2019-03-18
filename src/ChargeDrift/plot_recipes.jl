@recipe function f(t::Trajectory; showlabel = true, myscale = 1.0)
    legendfont --> 15

    if t.indv_charge>0
        c --> :green
        showlabel == true ? label --> "h+" : label --> ""
    else
        c --> :red
        showlabel == true ? label --> "e-" : label --> ""
    end
    lw := 3
    @series begin
        t.x_path * myscale, t.y_path * myscale, t.z_path * myscale
      end
    @series begin
        label := ""
        st:= :scatter
        [t.endposition[1]* myscale],[t.endposition[2]* myscale],[t.endposition[3]* myscale]
    end
end

@recipe function f(dp::DriftPath{T}; showlabel = true, scaling = 1.0) where T <: Real
    legendfont --> 15
    lw --> 3
    @series begin
        showlabel == true ? label --> "e-" : label --> ""
        c-->:red
        dp.e_path .*scaling
    end
    @series begin
        showlabel == true ? label --> "h+" : label --> ""
        c-->:green
        dp.h_path .* scaling
    end
    @series begin
        label := ""
        st:= :scatter
        c --> :red
        dp.e_path[end].* scaling
    end
    @series begin
        label := ""
        st:= :scatter
        c --> :green
        dp.h_path[end].* scaling
    end
end
@recipe function f(dps::AbstractVector{DriftPath{T}},scaling = 1.0) where T <:Real
    for i in eachindex(dps)
        @series begin
            scaling --> scaling
            i==1 ? showlabel --> true : showlabel --> false
            dps[i]
        end
    end
end
@recipe function f(p::Pulse)
    tickfont --> 20
    guidefont --> 20
    legendfont --> 10
    lw --> 4
    ylabel--> "Signal [keV]"
    xlabel --> "time [ns]"
    @series begin
        label := "e- contribution"
        c := :red
        p.electron_contribution
    end
    @series begin
        label := "h+ contribution"
        c := :green
        p.hole_contribution
    end
    @series begin
        label := "signal"
        c := :black
        p.signal
    end
end
@recipe function f(cde::ChargeDriftEvent, d::SolidStateDetector, show_pulses=true)

    if size(cde.signal,1)==1
        height = 300
        mywidths = [1/3,2/3]
        size --> (1200,500)
        layout --> (1, 2) #grid(1,2,widths=mywidths)
    else
        n_wps = size(cde.signal,1)
        width = 800
        length = width+n_wps*width/2
        size  -->  (width,length)
        myheights = [width/length]
        for i in 1:n_wps
            push!(myheights,(width/2)/length)
        end
        # layout  -->  (1, 2) #grid(n_wps+1,1,heights=myheights)
        layout --> (1+size(cde.signal,1),1)
    end
    @series begin
        aspect_ratio --> 1
        subplot --> 1
        # xlabel --> "x / mm"
        d
    end
    @series begin
        subplot --> 1
        scaling --> 1
        cde.drift_paths
    end

    for i in eachindex(cde.signal)
        subplot := i+1
        @series begin
            c --> :black
            lw --> 2
            label --> ""
            ylabel --> "Signal"
            xlabel --> "time / ns"
            [j for j in 1:size(cde.signal[1],1)], cde.signal[i]
        end
    end
end

@recipe function f(e::Event;coloring=[], labeling=[],show_pulses=true)
    if show_pulses==true
        n_wps = size(e.weighting_potentials,1)
        width = 800
        length = width+n_wps*width/2
        size  -->  (width,length)
        myheights = [width/length]
        for i in 1:n_wps
            push!(myheights,(width/2)/length)
        end

        layout  -->  (n_wps+1,1) #grid(n_wps+1,1,heights=myheights)
        @series begin
            subplot := 1
            coloring --> coloring
            labeling --> labeling
            e.detector
        end
        for itr in 1:e.n_sites
            @series begin
                itr == 1 ? showlabel --> true : showlabel --> false
                subplot := 1
                myscale := 1/e.detector.geometry_unit_factor
                e.trajectories_e[itr]
            end
            @series begin
                itr == 1 ? showlabel --> true : showlabel --> false
                subplot := 1
                myscale := 1/e.detector.geometry_unit_factor
                e.trajectories_h[itr]
            end
        end
        for iwp in 1:n_wps
            @series begin
                subplot := 1+iwp
                e.pulses[iwp]
            end
        end
    else
        width,length = 800, 800
        size  -->  (width,length)
        @series begin
            subplot := 1
            coloring --> coloring
            labeling --> labeling
            e.detector
        end

        for itr in 1:e.n_sites
            @series begin
                itr == 1 ? showlabel --> true : showlabel --> false
                subplot := 1
                myscale := 1/e.detector.geometry_unit_factor
                e.trajectories_e[itr]
            end
            @series begin
                itr == 1 ? showlabel --> true : showlabel --> false
                subplot := 1
                myscale := 1/e.detector.geometry_unit_factor
                e.trajectories_h[itr]
            end
        end
    end
end
