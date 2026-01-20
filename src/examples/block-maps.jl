include("../BlockMaps.jl")
using .BlockMaps

# To define block maps, use the `BlockMap` constructor.
"exchange between `1` and `2`"
flip_aut = BlockMap([2;;], [2;;], Dict(
    (1,) => 2,
    (2,) => 1,
), 0)

"shift map on the full 2-shift"
shift_aut = BlockMap([2;;], [2;;], Dict(
    (1,) => 1,
    (2,) => 2,
), -1)

# There is also a convenience function to define the shift maps.
shift([2;;])

# And a convenience function to define the compound marker endomorphisms on full 2-shift.
flip_aut = marker_endomorphism("*")

"composition of the shift and the flip"
shift_flip = flip_aut ∘ shift_aut

# Applying a block map to a sequence is done by calling it like a function.
shift_flip([1, 2, 1, 1, 2, 2])
# (2, 1, 2, 2, 1, 1)

# To calculate the dimension representation of a block map, as the first step use `_kichens_StrongShiftEquivalence`.
sse_flip, l_flip = BlockMaps._kichens_StrongShiftEquivalence(flip_aut)
display(sse_flip)
# StrongShiftEquivalence with lag 3:
#           ⎡1⎤        ⎡0 1⎤ ⎡1 1⎤        ⎡1⎤           
#     [1 1] ⎣1⎦  ⎡1 1⎤ ⎣1 0⎦ ⎣1 1⎦  ⎡1 1⎤ ⎣1⎦ [1 1]     
# [2]───────────▶︎⎣1 1⎦─────────────▶︎⎣1 1⎦───────────▶︎[2]

# Now we can compute the dimension representation.
dimension_representation(sse_flip, l_flip)
# 1×1 Matrix{Float64}:
#  1.0

"a example of non-inert automorphism"
scott_aut = BlockMap([2 1; 1 2], Dict(
    (1,) => 6,
    (2,) => 5,
    (3,) => 4,
    (4,) => 3,
    (5,) => 2,
    (6,) => 1,
), 0)

sse_scott, l_scott = BlockMaps._kichens_StrongShiftEquivalence(scott_aut)
dimension_representation(sse_scott, l_scott)
# 2×2 Matrix{Float64}:
#  0.0  1.0
#  1.0  2.22045e-16
