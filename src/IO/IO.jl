
function NamedTupleFromStruct(str)
    nt_struct = NamedTupleTools.ntfromstruct(str)
    nt_type = (type = typeof(str))
    return (nt_type = nt_type, nt_struct = nt_struct)
end

function StructFromNamedTuple(nt::NamedTuple)
    return NamedTupleTools.structfromnt(nt.nt_type, nt.nt_struct)
end

