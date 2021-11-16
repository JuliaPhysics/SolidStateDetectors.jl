@recipe function f(dp::EHDriftPath{T}; showlabel = true) where {T <: Real}
    legendfontsize --> 15
    linewidth --> 3
    seriestype --> :path3d
    @series begin
        showlabel ? label --> "e-" : label --> ""
        seriescolor --> :red
        [CartesianPoint{T}(path) for path in dp.e_path]
    end
    @series begin
        showlabel ? label --> "h+" : label --> ""
        seriescolor --> :green
        [CartesianPoint{T}(path) for path in dp.h_path]
    end
    @series begin
        label := ""
        seriestype := :scatter
        seriescolor --> :red
        CartesianPoint{T}(dp.e_path[end])
    end
    @series begin
        label := ""
        seriestype := :scatter
        seriescolor --> :green
        CartesianPoint{T}(dp.h_path[end])
    end
end
@recipe function f(dps::Vector{EHDriftPath{T}}) where {T <: Real}
    for i in eachindex(dps)
        @series begin
            showlabel --> i == 1
            dps[i]
        end
    end
end
