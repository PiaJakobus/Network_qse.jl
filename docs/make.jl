using Network_qse
using Documenter

push!(LOAD_PATH, "../src/")
makedocs(;
    modules=[Network_qse],
    authors="Pia Jakobus",
    repo="https://github.com/PiaJakobus/Network_qse.jl/blob/{commit}{path}#L{line}",
    sitename="Network QSE",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", nothing) == "true",
        canonical="https://PiaJakobus.github.io/Network_qse.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    branch="gh-pages",
    devbranch = "main",
    devurl = "stable",
    repo="github.com/PiaJakobus/Network_qse.jl.git",
)
