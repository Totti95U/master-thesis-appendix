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
        "Block maps" => "block-maps.md",
        "Cascade count" => "cascade-count.md",
        "Pruning" => "pruning.md",
        "Subshift lattice" => "subshift-lattice.md"
    ]
)