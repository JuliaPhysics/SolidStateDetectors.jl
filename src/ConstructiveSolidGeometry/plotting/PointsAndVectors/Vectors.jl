@recipe function f(::Type{Val{:vector}}, x, y, z; 
        headlength = 0.1, headwidth = 0.05 )
    origin = x
    vector = y
    target = origin + vector
    hl = headlength
    hw = headwidth
    T = eltype(origin)
    eX = CartesianVector{T}(1,0,0)
    eY = CartesianVector{T}(0,1,0)
    eZ = CartesianVector{T}(0,0,1)
    apxl = target - hl*vector - hw * normalize(vector × eX) * norm(vector) 
    apxr = target - hl*vector + hw * normalize(vector × eX) * norm(vector)
    apyl = target - hl*vector - hw * normalize(vector × eY) * norm(vector) 
    apyr = target - hl*vector + hw * normalize(vector × eY) * norm(vector)
    apzl = target - hl*vector - hw * normalize(vector × eZ) * norm(vector) 
    apzr = target - hl*vector + hw * normalize(vector × eZ) * norm(vector)

    label --> nothing
    linecolor --> :black
    xguide --> "x"
    # xunit --> internal_length_unit
    yguide --> "y"
    # yunit --> internal_length_unit
    zguide --> "z"
    # zunit --> internal_length_unit
    unitformat --> :slash
    if occursin("GRBackend", string(typeof(plotattributes[:plot_object].backend)))
        aspect_ratio --> 1.0
    end 
    x := [origin[1], target[1], apxl[1], target[1], apxr[1], target[1], apyl[1], target[1], apyr[1], target[1], apzl[1], target[1], apzr[1]]*internal_length_unit
    y := [origin[2], target[2], apxl[2], target[2], apxr[2], target[2], apyl[2], target[2], apyr[2], target[2], apzl[2], target[2], apzr[2]]*internal_length_unit
    z := [origin[3], target[3], apxl[3], target[3], apxr[3], target[3], apyl[3], target[3], apyr[3], target[3], apzl[3], target[3], apzr[3]]*internal_length_unit
    seriestype := :path3d
    ()
end

@recipe function f(v::CartesianVector{T}) where {T}
    @series begin
        seriestype --> :vector
        label --> nothing
        CartesianPoint{T}(0,0,0), v
    end
end