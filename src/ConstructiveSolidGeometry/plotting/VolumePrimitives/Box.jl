@recipe function f(b::Box)
    fs = faces(b)
    @series begin
        fs
    end
end
