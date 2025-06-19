using JLD2
using ScopedValues
using ParallelProcessingTools: read_files, write_files, CreateOrIgnore

const snapshot_var_path = ScopedValue("snapshots")

function load_or_generate(f_gen, varname::AbstractString, filename::AbstractString = "$varname.jld2")
    filepath = joinpath(snapshot_var_path[], filename)
    if isfile(filepath)
        @info "Loading $varname from \"$filepath\""
        return read_files(filepath, use_cache = false) do tmp_filepath
            JLD2.load(tmp_filepath, "obj")
        end
    else
        @info "Generating $varname and saving to \"$filepath\""
        obj = f_gen()

        write_files(
            filepath, mode = CreateOrIgnore(), use_cache = true,
            create_dirs = true
        ) do tmp_filepath
            JLD2.save(tmp_filepath, "obj", obj)
        end       

        return obj
    end
end


macro snapshot(expr)
    expr.head == :(=) || throw(ArgumentError("snapshot macro expects an assignment expression"))
    nonesc_varname = expr.args[1]
    nonesc_varname isa Symbol || throw(ArgumentError("snapshot macro expects a variable name on the left-hand side"))
    var = esc(nonesc_varname)
    varsym = QuoteNode(nonesc_varname)
    body = esc(expr.args[2])
    return quote
        varname = String($varsym)
        filename = joinpath(varname)
        $var = load_or_generate(() -> $body, filename)
    end
end
