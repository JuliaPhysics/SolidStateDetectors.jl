include("SigGenInterface.jl")
include("ParseConfigFiles.jl")

NamedTuple(::Missing) = (object = missing,)
NamedTuple(d::Dict) = (dict_json_string = json(d),)
Base.convert(T::Type{NamedTuple}, x::Dict) = T(x)
Dict(nt::NamedTuple) = JSON.parse(nt.dict_json_string)




