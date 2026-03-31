using Literate

# Config for the execution of notebooks
cfg = Dict(
    "credit" => false,
    "execute" => true,
    "continue_on_error" => false,
)

@info "Building to $(pwd())"
Literate.markdown("travail.jl", pwd(); config=cfg, flavor=Literate.CommonMarkFlavor())
