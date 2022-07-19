include("Mirror_Symmetry.jl")


#Converts all symmetry definitions inside nt into respective type
function Symmetry(nt::NamedTuple)
    out = NamedTuple()
    for sym_object in keys(nt)
        out[sym_object] = []
        for sym_type in keys(nt[key])
            if sym_type == "mirror_symmetry"
                append!(out[sym_object], MirrorSymmetry(nt[sym_object][sym_type]))
            end
        end
    end
    return out
end
 
