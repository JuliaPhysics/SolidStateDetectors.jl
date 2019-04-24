@recipe function f(dp::DriftPath{T}; showlabel = true, scaling = 1.0) where {T <: Real}
    legendfont --> 15
    lw --> 3
    @series begin
        showlabel == true ? label --> "e-" : label --> ""
        c --> :red
        [CartesianPoint{T}(path * scaling) for path in dp.e_path]
    end
    @series begin
        showlabel == true ? label --> "h+" : label --> ""
        c --> :green
        [CartesianPoint{T}(path * scaling) for path in dp.h_path]
    end
    @series begin
        label := ""
        st := :scatter
        c --> :red
        CartesianPoint{T}(dp.e_path[end] * scaling)
    end
    @series begin
        label := ""
        st := :scatter
        c --> :green
        CartesianPoint{T}(dp.h_path[end] * scaling)
    end
end
@recipe function f(dps::Vector{DriftPath{T}}, scaling = 1.0) where {T <: Real}
    for i in eachindex(dps)
        @series begin
            scaling --> scaling
            i == 1 ? showlabel --> true : showlabel --> false
            dps[i]
        end
    end
end
