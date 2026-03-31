using Literate

# Config for the execution of notebooks
cfg = Dict(
    "credit" => false,
    "execute" => true,
    "continue_on_error" => true,
)

# Replace the included files
function replace_includes(str)

    included = filter(endswith(".jl"), readdir("code", join=true))

    for ex in included
        content = read(ex, String)
        str = replace(str, "include(\"$(ex)\")" => content)
    end
    return str
end

@info "Building to $(pwd())"
Literate.markdown("travail.jl", pwd(); config=cfg, flavor=Literate.CommonMarkFlavor(), preprocess=replace_includes)
