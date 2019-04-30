
function NamedTupleFromStruct(str)
    nt_struct = NamedTupleTools.ntfromstruct(str)
    nt_type = (type = typeof(str))
    return (nt_type = nt_type, nt_struct = nt_struct)
end

function StructFromNamedTuple(nt::NamedTuple)
    return NamedTupleTools.structfromnt(nt.nt_type, nt.nt_struct)
end

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
Base.convert(T::Type{Dict}, x::NamedTuple) = T(x)




