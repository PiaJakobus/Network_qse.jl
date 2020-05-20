using Documenter
include("../src/Network_qse.jl")
using .Network_qse

makedocs(
    sitename = "Network_qse.jl",
    format = Documenter.HTML(),
    modules = [Network_qse]
)

# Documenter can also automatically deploy documentation to gh-pages.
# See "Hosting Documentation" and deploydocs() in the Documenter manual
# for more information.
deploydocs(
    repo = "github.com/PiaJakobus/Network_qse.jl.git", 
    devbranch = "master",
    devurl = "stable"
)
