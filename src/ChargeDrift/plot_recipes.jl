@recipe function f(dp::EHDriftPath{T}; showlabel = true, scaling = 1.0) where {T <: Real}
    legendfont --> 15
    linewidth --> 3
    @series begin
        showlabel == true ? label --> "e-" : label --> ""
        seriescolor --> :red
        [CartesianPoint{T}(path * scaling) for path in dp.e_path]
    end
    @series begin
        showlabel == true ? label --> "h+" : label --> ""
        seriescolor --> :green
        [CartesianPoint{T}(path * scaling) for path in dp.h_path]
    end
    @series begin
        label := ""
        seriestype := :scatter
        seriescolor --> :red
        CartesianPoint{T}(dp.e_path[end] * scaling)
    end
    @series begin
        label := ""
        seriestype := :scatter
        seriescolor --> :green
        CartesianPoint{T}(dp.h_path[end] * scaling)
    end
end
@recipe function f(dps::Vector{EHDriftPath{T}}, scaling = 1.0) where {T <: Real}
    for i in eachindex(dps)
        @series begin
            scaling --> scaling
            i == 1 ? showlabel --> true : showlabel --> false
            dps[i]
        end
    end
end
