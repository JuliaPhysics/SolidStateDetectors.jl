function NamedTuple(::Missing)
    return (object = missing,)
end


function NamedTuple(d::Dict) 
    return (dict_json_string = json(d),)
end
Base.convert(T::Type{NamedTuple}, x::Dict) = T(x)

function Dict(nt::NamedTuple)
    JSON.parse(nt.dict_json_string)
end
# Base.convert(T::Type{Dict}, x::NamedTuple) = T(x) # is now defined in NamedTupleTools-Package




