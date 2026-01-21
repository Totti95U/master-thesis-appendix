# Plot functions for subshift lattice
using GraphPlot

include("../subshift-lattice.jl")

"example useage of subshift relation ⊏"
[ [1, 0, 1, 0, 0], [1, 1, 1, 0, 0] ] ⊏ [ [0, 1, 0, 1, 0, 0], [0, 1, 1, 1, 0, 0] ]
# true

"example usage of subshift lattice plotting"
V = [
    [ [1, 0, 1, 0, 0], [1, 1, 1, 0, 0] ],
    [ [0, 1, 0, 1, 0, 0], [0, 1, 1, 1, 0, 0] ],
    [ [0, 0, 1, 0, 1, 0, 0], [0, 0, 1, 1, 1, 0, 0] ],
    [
        [0, 0, 1, 0], [0, 1, 1, 0],
        [1, 0, 0, 1, 1, 1], [1, 0, 1, 1, 1, 1],
    ],
    [ 
        [0, 0, 1, 0], [0, 1, 1, 0],
        [1, 1, 0, 0, 1, 1, 1, 0], [1, 1, 0, 1, 1, 1, 1, 0],
    ],
    [ [0, 0, 1, 0], [0, 1, 1, 0] ],
    [
        [1, 0, 0, 1, 0], [1, 0, 1, 1, 0],
        [0, 0, 0, 0, 1, 0], [0, 0, 0, 1, 1, 0],
        [1, 0, 0, 1, 1, 1], [1, 0, 1, 1, 1, 1],
    ],
    [ [1, 0, 0, 1, 0], [1, 0, 1, 1, 0] ],
    Word[] #full shift (`Word` is an alias for `Vector{Int}`)
]

shift_hasse_diagram(V)