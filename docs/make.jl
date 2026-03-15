using Documenter

include("../src/BlockMaps.jl")
include("../src/cascade-count.jl")
include("../src/Pruning.jl")
include("../src/subshift-lattice.jl")
makedocs(
    sitename = "Appendix of my master thesis", 
    remotes=nothing,
    pages = Any[
        "Home" => "index.md",
        "Guide" => "guide.md",
        "API" => [
            "BlockMaps.jl" => "block-maps-api.md",
            "cascade-count.jl" => "cascade-count-api.md",
            "Pruning.jl" => "pruning-api.md",
            "subshift-lattice.jl" => "subshift-lattice-api.md",
        ]
    ]
)