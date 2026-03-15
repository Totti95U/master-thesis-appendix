# Implement algorithm for finding primary pruning fronts of given real Hénon map
# See: Shudo & Hagiwara (2004), "An algorithm to prune the area-preserving Hénon map" from Journal of Physics. A. Mathematical and General

using CairoMakie, LinearAlgebra, SpatialIndexing

"working accuracy alias:"
const MyFloat = Float64

"Custom type for 2D points"
const Point2o = Point{2, MyFloat}

# make Point2o work with isfinite
Base.isfinite(pt::Point2o) = all(isfinite, pt)

"""
    HenonMap(a, b)

Create a Hénon map with parameters `a` and `b` defined by

```math
H_{a, b}(x, y) = (-x^2 + b y + a, x)
```

The map is area-preserving if ``|b| = 1``.
When ``b = -1``, the map is orientation-preserving.
"""
struct HenonMap
    a::MyFloat
    b::MyFloat
end

function forward(hm::HenonMap, x::Point2o)
    new_x = -x[1]^2 + hm.b * x[2] + hm.a
    new_y = x[1]
    return Point2o(new_x, new_y)
end

function backward(hm::HenonMap, x::Point2o)
    new_x = x[2]
    new_y = (x[1] + x[2]^2 - hm.a) / hm.b
    return Point2o(new_x, new_y)
end

"""
    jacobian(hm::HenonMap, x::Point2o)

Return the Jacobian matrix of the Hénon map `hm` at point `x`.
"""
function jacobian(hm::HenonMap, x::Point2o)
    return [
        -2 * x[1]  hm.b;
            1        0
    ]
end

"""
    eigen(hm::HenonMap, pt::Point2o) -> (eigenvalues::Vector{MyFloat}, eigenvectors::Matrix{MyFloat})

Return the eigenvalues and eigenvectors of the Jacobian of the Hénon map `hm` at point `pt`.
"""
function LinearAlgebra.eigen(hm::HenonMap, pt::Point2o)
    x = pt[1]
    d = sqrt(hm.b + x^2)
    λs = [-x + d, -x - d]
    vs = [Point2o(λ, 1) for λ in λs]
    return λs, hcat(vs...)
end

"""
    hyperbolic_fixed_points(hm::HenonMap)

Return the hyperbolic fixed points of the Hénon map `hm`.
"""
function hyperbolic_fixed_points(hm::HenonMap)
    disc = (hm.b - 1)^2 + 4 * hm.a
    if disc < 0
        return nothing
    end
    sqrt_disc = sqrt(disc)
    x1 = ((hm.b - 1) + sqrt_disc) / 2
    x2 = ((hm.b - 1) - sqrt_disc) / 2
    p1 = Point2o(x1, x1)
    p2 = Point2o(x2, x2)
    fixpts = [p1, p2]
    hyperbolic_pts = []
    for p in fixpts
        J = jacobian(hm, p)
        specs = abs.(eigvals(J))
        if (specs[1] > 1 && specs[2] < 1) || (specs[1] < 1 && specs[2] > 1)
            push!(hyperbolic_pts, p)
        end
    end

    return hyperbolic_pts
end

"""
    image_of(f::Function, xs::Vector{Point2o}; tolerance)

Return the image of a set of points `xs` under function `f`, adaptively refining `xs` to ensure
that the image is continuous within `tolerance`.
"""
function image_of(f::Function, xs::Vector{Point2o}; tolerance)
    is_image_discontinued = true

    norm(x::Point2o) = sqrt(x[1]^2 + x[2]^2)
    
    while is_image_discontinued
        is_image_discontinued = false

        num_new_midpoints = 0
        for_loop_num = length(xs) - 1

        for i in 1:for_loop_num
            ii = i + num_new_midpoints
            if norm(f(xs[ii]) - f(xs[ii+1])) > tolerance
                is_image_discontinued = true
                # insert midpoint between xs[i] and xs[i+1] at position i+1 in tmp_xs
                midpoint = (xs[ii] + xs[ii+1]) / 2
                insert!(xs, ii+1, midpoint)
                num_new_midpoints += 1
            end
        end
    end

    return [f(x) for x in xs]
end

"""
    manifolds(hm::HenonMap; num_iterations=16, delta=1e-5, tolerance=1e-3, fixed_pt::Point2o)

Return the stable and unstable manifolds of the hyperbolic fixed points of the Hénon map `hm`.
If `fixed_pt_id` is provided, only compute manifolds for the specified fixed points (1-based index).
"""
function manifolds(hm::HenonMap; num_iterations=16, delta::MyFloat=MyFloat(1e-5), tolerance=1e-3, fixed_pt::Point2o)
    stable_manifolds = Vector{Point2o}[]
    unstable_manifolds = Vector{Point2o}[]

    eigenvals, eigenvecs = eigen(hm, fixed_pt)
    for (i, λ) in enumerate(eigenvals)
        branches = Vector{Point2o}[] # the set of segments of the manifold
        xs = [Point2o(0, 0), delta * Point2o(eigenvecs[:, i])] .+ fixed_pt # initial segment of the manifold
        tmp_f = abs(λ) < 1 ? x -> backward(hm, x) : x -> forward(hm, x)
        push!(branches, xs)
        for _ in 1:num_iterations
            old_branches = copy(branches)
            branches = Vector{Point2o}[]
            for xs in old_branches
                xs = image_of(tmp_f, xs; tolerance=tolerance)
                # trim the manifold to only include points within [-5, 5] x [-5, 5]
                xs = filter(x -> all(-5 .<= x .<= 5), xs)
                # divide xs at the gaps larger than tolerance
                cut_idx = Int[1]
                for i in axes(xs[1:end-1], 1)
                    if norm(xs[i+1] - xs[i]) > 10 * tolerance
                        push!(cut_idx, i+1)
                    end
                end
                if isempty(cut_idx)
                    push!(branches, xs)
                    continue
                end
                for i in axes(cut_idx[1:end-1], 1)
                    push!(branches, xs[cut_idx[i]:cut_idx[i+1]-1])
                end
                push!(branches, xs[cut_idx[end]:end])
            end
        end

        if abs(λ) < 1
            append!(stable_manifolds, branches)
        else
            append!(unstable_manifolds, branches)
        end
    end

    return stable_manifolds, unstable_manifolds
end

"""
    primary_branch(hm::HenonMap, partition::Function; delta=1e-5, tolerance=1e-3, fixed_pt::Point2o)

Return the primary branches of the stable and unstable manifolds of the Hénon map `hm` at `fixed_pt`.
The primary branch is defined as the branch that remains in the same partition as the fixed point under iteration.
"""
function primary_branch(hm::HenonMap, partition::Function; delta::MyFloat=MyFloat(1e-5), tolerance=1e-3, fixed_pt::Point2o)
    eigenvals, eigenvecs = eigen(hm, fixed_pt)
    stable_branch = Point2o[]
    unstable_branch = Point2o[]
    for (i, λ) in enumerate(eigenvals)
        xs = [Point2o(0, 0), delta * Point2o(eigenvecs[:, i])] .+ fixed_pt # initial segment of the manifold
        tmp_f = abs(λ) < 1 ? x -> backward(hm, x) : x -> forward(hm, x)
        for _ in 1:16
            xs = image_of(tmp_f, xs; tolerance=tolerance)
            for (i, x) in enumerate(xs)
                if partition(x) != partition(fixed_pt)
                    xs = xs[1:i-1]
                    break
                end
            end
        end
        if abs(λ) < 1
            append!(stable_branch, xs)
        else
            append!(unstable_branch, xs)
        end
    end
    return stable_branch, unstable_branch
end

"""
    _seg_intersection(p1, p2, q1, q2) -> Union{Point2o, Nothing}

Check intersection between segment p1-p2 and q1-q2.
Return intersection point if it exists, otherwise `nothing`.
"""
function _seg_intersection(p1::Point2o, p2::Point2o, q1::Point2o, q2::Point2o)
    x1, y1 = p1
    x2, y2 = p2
    x3, y3 = q1
    x4, y4 = q2

    # Solve determinant
    denom = (x1-x2)*(y3-y4) - (y1-y2)*(x3-x4)
    if iszero(denom)
        return nothing # parallel or collinear
    end
    px = ((x1*y2 - y1*x2)*(x3-x4) - (x1-x2)*(x3*y4 - y3*x4)) / denom
    py = ((x1*y2 - y1*x2)*(y3-y4) - (y1-y2)*(x3*y4 - y3*x4)) / denom

    # TODO: 端点での交差も含めた方がいいかもしれない
    if min(x1,x2) - 1e-9 <= px <= max(x1,x2) + 1e-9 &&
       min(y1,y2) - 1e-9 <= py <= max(y1,y2) + 1e-9 &&
       min(x3,x4) - 1e-9 <= px <= max(x3,x4) + 1e-9 &&
       min(y3,y4) - 1e-9 <= py <= max(y3,y4) + 1e-9
        return Point2o(px, py)
    end
    return nothing
end

"""
    _poly_segments(lines::Vector{Vector{Point2o}}) -> Vector{Tuple{Point2o,Point2o,Int,Int}}

Return all line segments from given polylines.
Each segment is (p1, p2, poly_idx, seg_idx).
"""
function _poly_segments(lines::Vector{Vector{Point2o}})
    segs = Vector{Tuple{Point2o,Point2o,Int,Int}}()
    for (pi, poly) in enumerate(lines)
        for si in 1:length(poly)-1
            push!(segs, (poly[si], poly[si+1], pi, si))
        end
    end
    segs
end

"""
    _item_index(item) -> Int

Extract integer index out of whatever `intersects_with` yields.
Handles: direct integer, wrapper with `.id`, wrapper with `.value`, or `id(item)`.
"""
function _item_index(item)
    # direct integer
    if item isa Integer
        return Int(item)
    end

    # try common accessors; different versions return different shapes
    try
        return Int(getfield(item, :id))
    catch
    end
    try
        return Int(getfield(item, :value))
    catch
    end
    try
        return Int(id(item))
    catch
    end

    error("cannot extract stored id from spatial index item: $(typeof(item))")
end

"""
    _build_index(segs)

Build an RTree storing (key, value) = (segment_index, segment_index).
Use SpatialIndexing.Rect/Vec (not GeometryBasics types).
"""
function _build_index(segs::Vector{Tuple{Point2o,Point2o,Int,Int}})
    # construct empty tree: specify numeric coordinate type and (KeyType, ValueType)
    # here we store Int keys and Int values
    tree = RTree{MyFloat,2}(Int, Int)

    for (i, (p1, p2, _, _)) in enumerate(segs)
        x1, y1 = MyFloat(p1[1]), MyFloat(p1[2])
        x2, y2 = MyFloat(p2[1]), MyFloat(p2[2])
        xmin, xmax = min(x1,x2), max(x1,x2)
        ymin, ymax = min(y1,y2), max(y1,y2)
        # create SpatialIndexing rect (NOT GeometryBasics HyperRectangle)
        rect = SpatialIndexing.Rect((xmin, ymin), (xmax, ymax))
        # insert with (key, value). We choose to store the segment's index as both key & value.
        insert!(tree, rect, i, i)
    end

    return tree
end

"""
    intersections(lines1, lines2) -> Vector{Tuple{Int,Int,Int,Int,Point2o}}

Compute all intersections between polylines in `lines1` and `lines2`.
Returns tuples: (poly1_idx, seg1_idx, poly2_idx, seg2_idx, intersection_point).
"""
function intersections(lines1::Vector{Vector{Point2o}}, lines2::Vector{Vector{Point2o}})
    segs1 = _poly_segments(lines1)
    segs2 = _poly_segments(lines2)
    tree2 = _build_index(segs2)

    result = Vector{Tuple{Int,Int,Int,Int,Point2o}}()

    for (p1,p2,pi1,si1) in segs1
        x1, y1 = MyFloat(p1[1]), MyFloat(p1[2])
        x2, y2 = MyFloat(p2[1]), MyFloat(p2[2])
        xmin, xmax = min(x1,x2), max(x1,x2)
        ymin, ymax = min(y1,y2), max(y1,y2)
        qrect = SpatialIndexing.Rect((xmin, ymin), (xmax, ymax))

        for item in intersects_with(tree2, qrect)
            idx = _item_index(item)            # index into segs2
            q1,q2,pi2,si2 = segs2[idx]
            ip = _seg_intersection(p1, p2, q1, q2)
            if ip !== nothing
                push!(result, (pi1, si1, pi2, si2, ip))
            end
        end
    end
    return result
end

"""
    HomoclinicCode(backward::Vector{Int}, forward::Vector{Int})

A type representing the symbolic code of a homoclinic point, with forward and backward components.

# Example
```jldoctest pruning
julia> HomoclinicCode([1, 0], [1,0,1])
5-element HomoclinicCode:
(0ᵒᵒ)10.101(0ᵒᵒ)

julia> HomoclinicCode("10.101") # string constructor
5-element HomoclinicCode:
(0ᵒᵒ)10.101(0ᵒᵒ)
```
"""
struct HomoclinicCode
    backward::Vector{Int}
    forward::Vector{Int}

    function HomoclinicCode(backward::Vector{Int}, forward::Vector{Int})
        fwd_end = findlast(!=(0), forward)
        bwd_start = findfirst(!=(0), backward)

        new_fwd = fwd_end === nothing ? Int[] : forward[1:fwd_end]
        new_bwd = bwd_start === nothing ? Int[] : backward[bwd_start:end]
        new(new_bwd, new_fwd)
    end
end

"Constructor from a tuple of two vectors. The first element will be the backward code, the second will be the forward code."
HomoclinicCode(c::Tuple{Vector{Int}, Vector{Int}}) = HomoclinicCode(c[1], c[2])

"""
    HomoclinicCode(s::String)

Construct a `HomoclinicCode` from a string representation like "1000.0001".
The string should contain exactly one '.' character separating the backward and forward parts.
"""
function HomoclinicCode(s::String)
    if !occursin('.', s)
        error("invalid code string: missing '.' separator")
    end
    parts = split(s, '.')
    if length(parts) != 2
        error("invalid code string: should contain exactly one '.' separator")
    end
    bwd_str, fwd_str = parts
    bwd_vec = [parse(Int, c) for c in bwd_str]
    fwd_vec = [parse(Int, c) for c in fwd_str]
    return HomoclinicCode(bwd_vec, fwd_vec)
end

function Base.:(==)(hc1::HomoclinicCode, hc2::HomoclinicCode)
    return hc1.backward == hc2.backward && hc1.forward == hc2.forward
end

Base.isequal(hc1::HomoclinicCode, hc2::HomoclinicCode) = hc1 == hc2

Base.hash(hc::HomoclinicCode, h::UInt) = hash((hc.backward, hc.forward), h)

function Base.show(io::IO, hc::HomoclinicCode)
    print(io, "HomoclinicCode(")
    print(io, join(hc.backward, ""))
    print(io, ".")
    print(io, join(hc.forward, ""))
    print(io, ")")
end

function Base.show(io::IO, ::MIME"text/plain", hc::HomoclinicCode)
    println(io, "$(length(hc))-element HomoclinicCode:")
    print(io, "(0ᵒᵒ)")
    print(io, join(hc.backward, ""))
    print(io, ".")
    print(io, join(hc.forward, ""))
    print(io, "(0ᵒᵒ)")
end

Base.length(hc::HomoclinicCode) = length(hc.backward) + length(hc.forward)

"""
    Base.getindex(hc::HomoclinicCode, i::Int)
"""
function Base.getindex(hc::HomoclinicCode, i::Int) 
    if i > length(hc.forward)-1 || i < -length(hc.backward)
        return 0
    end

    return i < 0 ? hc.backward[end+1+i] : hc.forward[i+1]
end
Base.getindex(hc::HomoclinicCode, I::AbstractVector{<:Integer}) = [hc[i] for i in I]
Base.getindex(hc::HomoclinicCode, r::UnitRange{<:Integer}) = [hc[i] for i in r]

"""
    ≺(s::HomoclinicCode, t::HomoclinicCode) -> Tuple{Bool, Bool}

Unimodal order on symbolic codes.
"""
function ≺(s::HomoclinicCode, t::HomoclinicCode)
    bwd_res = false
    if s.backward == t.backward
        bwd_res = false
    else
        unequi_idx = 0
        for i in -1:-1:min(-length(s.backward), -length(t.backward))
            if s[i] != t[i]
                unequi_idx = i
                break
            end
        end
        if iseven(sum(s[unequi_idx+1:-1]))
            bwd_res = s[unequi_idx] < t[unequi_idx]
        else
            bwd_res = s[unequi_idx] > t[unequi_idx]
        end
    end
    fwd_res = false
    if s.forward == t.forward
        fwd_res = false
    else
        unequi_idx = 0
        for i in 0:max(length(s.forward)-1, length(t.forward)-1)
            if s[i] != t[i]
                unequi_idx = i
                break
            end
        end
        if iseven(sum(s[0:unequi_idx-1]))
            fwd_res = s[unequi_idx] < t[unequi_idx]
        else
            fwd_res = s[unequi_idx] > t[unequi_idx]
        end
    end
    return (bwd_res, fwd_res)
end

"""
    trim(hc::HomoclinicCode) -> HomoclinicCode

Trim leading and trailing zeros from the symbolic code.

So, the returned code has a form `1 s₋ₕ ... s₋₁ . s₀ s₁ ... sₖ 1`.
"""
function trim(hc::HomoclinicCode)
    fwd_end = findlast(!=(0), hc.forward)
    bwd_start = findfirst(!=(0), hc.backward)

    new_fwd = fwd_start === nothing ? Int[] : hc.forward[1:fwd_end]
    new_bwd = bwd_start === nothing ? Int[] : hc.backward[bwd_start:end]

    return HomoclinicCode(new_fwd, new_bwd)
end

"""
Convert a binary code to its lexicographically equivalent form. This is needed for correct ordering in the symbolic plane. This is the inverse of `_lexi_to_knead`.
The first digit is the least significant bit in the lexicographical order, and the last digit is the most significant bit.
"""
function _knead_to_lexi(code::Vector{Int})
    w = copy(code)
    acc = accumulate(+, w)
    for i in axes(w[1:end-1], 1)
        w[i+1] = acc[i] % 2 == 0 ? w[i+1] : 1 - w[i+1]
    end
    return reverse(w)
end

"Convert a binary code from lexicographical form back to kneading form. This is the inverse of `_knead_to_lexi`."
function _lexi_to_knead(code::Vector{Int})
    tmp_code = reverse(copy(code))
    w = similar(code)
    # acc = accumulate(+, code)
    w[1] = tmp_code[1]
    for i in axes(w[1:end-1], 1)
        w[i+1] = sum(w[1:i]) % 2 == 0 ? tmp_code[i+1] : 1 - tmp_code[i+1]
    end
    return w
end

"squared Euclidean distance between two 2D points/tuples"
@inline function sqdist(a::Point2o, b::Point2o)
    dx = MyFloat(a[1]) - MyFloat(b[1])
    dy = MyFloat(a[2]) - MyFloat(b[2])
    return dx*dx + dy*dy
end

"""
squared distance from point p to segment a-b and the closest point on the segment
returns (d2, cx, cy, t) where t in [0,1] is the barycentric coordinate along a->b
"""
@inline function seg_point_sqdist(a::Point2o, b::Point2o, p::Point2o)
    ax = MyFloat(a[1]); ay = MyFloat(a[2])
    bx = MyFloat(b[1]); by = MyFloat(b[2])
    px = MyFloat(p[1]); py = MyFloat(p[2])

    vx = bx - ax
    vy = by - ay
    wx = px - ax
    wy = py - ay

    denom = vx*vx + vy*vy
    if denom == 0.0
        # degenerate segment: a == b
        cx, cy = ax, ay
        dx = px - cx
        dy = py - cy
        return (dx*dx + dy*dy, cx, cy, 0.0)
    end

    t = (wx*vx + wy*vy) / denom
    if t <= 0.0
        cx, cy = ax, ay
        dx = px - cx; dy = py - cy
        return (dx*dx + dy*dy, cx, cy, 0.0)
    elseif t >= 1.0
        cx, cy = bx, by
        dx = px - cx; dy = py - cy
        return (dx*dx + dy*dy, cx, cy, 1.0)
    else
        cx = ax + t*vx
        cy = ay + t*vy
        dx = px - cx; dy = py - cy
        return (dx*dx + dy*dy, cx, cy, t)
    end
end

"""
closest_point_bruteforce(p::Point2o, poly::Vector{Point2o})

Returns (min_dist, closest_point::Point2o, seg_index::Int, t::Float64)
- min_dist is the Euclidean distance (not squared)
- seg_index is the 1-based index of the segment start: segment = (poly[seg_index], poly[seg_index+1])
- t is barycentric coordinate in [0,1] along that segment
"""
function closest_point_bruteforce(p::Point2o, poly::Vector{Point2o})
    n = length(poly)
    if n == 0
        error("polyline has no points")
    elseif n == 1
        return (sqrt(sqdist(poly[1], p)), Point2o(poly[1]), 0, 0.0)
    end

    best_d2 = typemax(MyFloat)
    best_pt = Point2o(0.0, 0.0)
    best_seg = 0
    best_t = 0.0

    @inbounds for i in 1:(n-1)
        a = poly[i]; b = poly[i+1]
        d2, cx, cy, t = seg_point_sqdist(a, b, p)
        if d2 < best_d2
            best_d2 = d2
            best_pt = Point2o(MyFloat(cx), MyFloat(cy))
            best_seg = i
            best_t = t
        end
    end

    return (sqrt(best_d2), best_pt, best_seg, best_t)
end

"""
    prepare_polyline(polyline::Vector{Point2o}) -> (segs, tree)

Prepare a polyline for batch closest point queries.
"""
function prepare_polyline(polyline::Vector{Point2o})
    # build segments and index
    segs = Vector{Tuple{Point2o,Point2o,Int,Int}}()
    for i in 1:(length(polyline)-1)
        push!(segs, (polyline[i], polyline[i+1], 1, i))
    end
    tree = _build_index(segs)
    return segs, tree
end

"""
    closest_point(pt::Point2o, seg, tree)

For point `pt`, find the closest point on the given `polyline`.

Arguments:
- `pt`: a query point (Point2o)
- segs: segments of the polyline, as returned by `prepare_polyline`
- tree: spatial index of the segments, as returned by `prepare_polyline`

Returns a tuple:
    `(dist::Float64, closest_point::Point2o, seg_idx::Int, t::Float64)`

where:
- `dist` is Euclidean distance
- `closest_point` is projection point on the polyline
- `seg_idx` is the starting index of the segment in `polyline` (segment = (polyline[seg_idx], polyline[seg_idx+1]))
- `t ∈ [0,1]` is barycentric coordinate along that segment.
"""
function closest_point(p::Point2o, segs, tree)
    px, py = MyFloat(p[1]), MyFloat(p[2])

    best_d2 = typemax(MyFloat)
    best_pt = Point2o(0,0)
    best_seg = 0
    best_t = 0.0

    if !isfinite(px) || !isfinite(py)
        @warn "query point has NaN or Inf coordinate"
        return (NaN, Point2o(NaN, NaN), NaN, NaN)
    end

    radius = 10
    while best_seg == 0 && radius < 1e6 # set maximum radius to avoid infinite loop
        qrect = SpatialIndexing.Rect((MyFloat(px-radius), MyFloat(py-radius)), (MyFloat(px+radius), MyFloat(py+radius)))
        for item in intersects_with(tree, qrect)
            idx = _item_index(item)
            a,b,_,si = segs[idx]
            d2,cx,cy,t = seg_point_sqdist(a,b,p)
            if d2 < best_d2
                best_d2 = d2
                best_pt = Point2o(cx,cy)
                best_seg = si
                best_t = t
            end
        end
        radius *= 2
    end

    return (sqrt(best_d2), best_pt, best_seg, best_t)
end

"""
    symbolic_encoding(hm::HenonMap, pts::Vector{Point2o}, partition::Function, psb::Vector{Point2o}, pub::Vector{Point2o}; iter=15) -> HomoclinicCode

Compute symbolic encoding of each point in `pts` under Hénon map `hm` using `partition` function.
The length of backward and forward codes are continued until iteration of the point hit the primary (un-)stable branch, `psb` (`pub`).
`partition` should map a `Point2o` to an integer symbol.
"""
function symbolic_encoding(hm::HenonMap, pts::Vector{Point2o}, partition::Function, psb::Vector{Point2o}, pub::Vector{Point2o}; iter=15)
    fwd_codes = Vector{Vector{Int}}(undef, length(pts))
    bwd_codes = Vector{Vector{Int}}(undef, length(pts))
    psb_segs, psb_tree = prepare_polyline(psb)
    pub_segs, pub_tree = prepare_polyline(pub)
    Threads.@threads for i in axes(pts, 1)
        pt = pts[i]
        fwd_code = Int[]
        x = pt
        for _ in 1:iter
            push!(fwd_code, partition(x))
            x = forward(hm, x)
            !isfinite(x) && break
            closest_point(x, psb_segs, psb_tree)[1] < 1e-3 && break
        end
        bwd_code = Int[]
        x = pt
        for _ in 1:iter
            x = backward(hm, x)
            !isfinite(x) && break
            closest_point(x, pub_segs, pub_tree)[1] < 1e-3 && break
            pushfirst!(bwd_code, partition(x))
        end

        fwd_codes[i] = fwd_code
        bwd_codes[i] = bwd_code
    end
    return [HomoclinicCode(bc, fc) for (bc, fc) in zip(bwd_codes, fwd_codes)]
end

"""
    in_symbolic_plane(symb_code::HomoclinicCode; ϵ=1e-2) -> Point2o

Transform a homoclinic code into a points in 2D symbolic plane [0, 1] x [0, 1].
"""
function in_symbolic_plane(symb_code::HomoclinicCode; ϵ=1e-2)
    bwd_code, fwd_code = copy(symb_code.backward), copy(symb_code.forward)
    pushfirst!(bwd_code, 0 )
    push!(fwd_code, 0)
    base = 2
    
    bwd_lexi = _knead_to_lexi(bwd_code |> reverse)
    fwd_lexi = _knead_to_lexi(fwd_code)
    y = sum(bwd_lexi[i] / base.^i for i in axes(bwd_lexi, 1))
    x = sum(fwd_lexi[i] / base.^i for i in axes(fwd_lexi, 1))
    y += sum(bwd_code) % 2 == 0 ? 0.0 : 1/base^length(bwd_code) - ϵ
    x += sum(fwd_code) % 2 == 0 ? 0.0 : 1/base^length(fwd_code) - ϵ
    return Point2o(x, y)
end

"""
    in_symbolic_plane_finite_tail(symb_code::HomoclinicCode, max_backward::Int, max_forward::Int) -> Point2o

Transform a homoclinic code into a points in 2D symbolic plane [0, 1] x [0, 1].
This version only considers up to `max_forward` symbols in the forward direction and `max_backward` symbols in the backward direction, treating the rest as zeros.
"""
function in_symbolic_plane_finite_tail(symb_code::HomoclinicCode, max_backward::Int, max_forward::Int)
    bwd_code, fwd_code = copy(symb_code[-max_backward:-1]), copy(symb_code[0:max_forward-1])
    base = 2
    
    bwd_lexi = _knead_to_lexi(bwd_code |> reverse)
    fwd_lexi = _knead_to_lexi(fwd_code)
    y = sum(bwd_lexi[i] / base.^i for i in axes(bwd_lexi, 1))
    x = sum(fwd_lexi[i] / base.^i for i in axes(fwd_lexi, 1))
    # ϵ = 1e-2
    # y += sum(bwd_code) % 2 == 0 ? 0.0 : 1/base^length(bwd_code) - ϵ
    # x += sum(fwd_code) % 2 == 0 ? 0.0 : 1/base^length(fwd_code) - ϵ
    return Point2o(x, y)
end

"""
    hat(hc::HomoclinicCode) -> HomoclinicCode

Return the symbolic code whose 0th symbol is flipped.
"""
function hat(hc::HomoclinicCode)
    if isempty(hc.forward)
        new_fwd = [1]
    else
        new_fwd = copy(hc.forward)
        new_fwd[1] = 1 - new_fwd[1]
    end
    return HomoclinicCode(hc.backward, new_fwd)
end

"""
    Nu(hc::HomoclinicCode) -> HomoclinicCode

Return the symbolic code whose secondary head symbol is flipped.
"""
function Nu(hc::HomoclinicCode)
    if length(hc.forward) < 2
        error("cannot apply Nu to a code with less than 2 forward symbols")
    end
    new_fwd = copy(hc.forward)
    new_fwd[end-1] = 1 - new_fwd[end-1]
    return HomoclinicCode(hc.backward, new_fwd)
end

"""
    Ns(hc::HomoclinicCode) -> HomoclinicCode

Return the symbolic code whose secondary tail symbol is flipped.
"""
function Ns(hc::HomoclinicCode)
    if length(hc.backward) < 2
        error("cannot apply Ns to a code with less than 2 backward symbols")
    end
    new_bwd = copy(hc.backward)
    new_bwd[2] = 1 - new_bwd[2]
    return HomoclinicCode(new_bwd, hc.forward)
end

"""
    Pu(h::Int, hc::HomoclinicCode) -> HomoclinicCode

Return the symbolic code which is pruned pair of `hc` by `u`-pruning at depth `h`.
"""
function Pu(h::Int, hc::HomoclinicCode)
    T = length(hc.forward)-1
    if h < T
        error("cannot apply Pu with h < length of forward code")
    end

    new_fwd = [copy(hc.forward)..., zeros(Int, h - T)...]
    new_fwd[end] = 1 - new_fwd[end]
    return HomoclinicCode(hc.backward, new_fwd)
end

"""
    Ps(t::Int, hc::HomoclinicCode) -> HomoclinicCode

Return the symbolic code which is pruned pair of `hc` by `s`-pruning at depth `s`.
"""
function Ps(t::Int, hc::HomoclinicCode)
    T = length(hc.backward)-1
    if t < T
        error("cannot apply Ps with t < length of backward code")
    end

    new_bwd = [zeros(Int, t - T)..., copy(hc.backward)...]
    new_bwd[1] = 1 - new_bwd[1]
    return HomoclinicCode(new_bwd, hc.forward)
end

"""
        primary_pruning_front(symb_codes::Vector{HomoclinicCode}) -> Vector{NTuple{4, HomoclinicCode}}

Given a list of homoclinic codes, compute the primary pruning front as a list of blocks.
Each block is represented as a tuple of four homoclinic codes `(a, b, c, d)`, where:
- `a` is the top-right corner of the block
- `b` is the bottom-right corner of the block
- `c` is the bottom-left corner of the block
- `d` is the top-left corner of the block
"""
function primary_pruning_front(symb_codes::Vector{HomoclinicCode})
    # 1. Setting Initial conditions
    a = HomoclinicCode([1,], [0, 1])
    j = 1
    t = 0
    b = HomoclinicCode(Int[], [0, 1])
    h = 2

    blocks = Vector{NTuple{4, HomoclinicCode}}() # (a, b, c, d) for each block D

    while true
        # 2. Find the depth of the block D
        top = a
        while true
            while true
                if (Ps(t, top) ∈ symb_codes) || (hat(Ps(t, top)) ∈ symb_codes)
                    t += 1
                else
                    break
                end
            end
            bot = Ps(t, top)
            top = Ns(bot)
            if (top ∈ symb_codes) ⊻ (top ≺ b)[1]
                b = bot
                break
            end
        end
        # 3. Find the width of the block D
        c = HomoclinicCode(Int[], Int[])
        while true
            if !((Pu(h, b) ∉ symb_codes) ⊻ (Pu(h, (Nu∘Pu)(h, b)) ∉ symb_codes))
                h += 1
                continue
            elseif Pu(h, b) ∉ symb_codes
                c = Pu(h, b)
                break
            else
                c = Pu(h, (Nu∘Pu)(h, b))
                break
            end
        end

        # 4. Complete the block D
        d = HomoclinicCode([1], c.forward)
        push!(blocks, (a, b, c, d))
        push!(blocks, (hat(a), hat(b), hat(c), hat(d)))

        # 5. setup for next iteration
        a = Nu(d)
        if a ∉ symb_codes && hat(a) ∉ symb_codes
            j += 1
        else
            break
        end
    end

    # 6. return the pruning front
    return blocks
end

"""
    forbidden_words(blocks::Vector{NTuple{4, HomoclinicCode}}) -> Vector{String}

Given a list of blocks representing the primary pruning front, compute the corresponding forbidden words.

# Example
```jldoctest pruning
julia> block = (HomoclinicCode("1.01"), HomoclinicCode("1001.01"), HomoclinicCode("1001.010001"), HomoclinicCode("1.010001"))
(HomoclinicCode(1.01), HomoclinicCode(1001.01), HomoclinicCode(1001.010001), HomoclinicCode(1.010001))

julia> forbidden_words([block])
1-element Vector{String}:
 "00101000"
```
"""
function forbidden_words(blocks::Vector{NTuple{4, HomoclinicCode}})
    words = Vector{String}()
    
    lexi_to_num(c::Vector{Int}) = sum(c[i] * 2^(i-1) for i in axes(c, 1))
    num_to_lexi(n::Int, len::Int) = [n >> (i-1) & 1 for i in 1:len]
    
    for (a, b, c, _) in blocks
        # compute backward part
        bwd_len = max(length(a.backward), length(b.backward))
        a_bwd_lexi = _knead_to_lexi(a[-1:-1:-bwd_len])
        b_bwd_lexi = _knead_to_lexi(b[-1:-1:-bwd_len])
        num_max = max(lexi_to_num(a_bwd_lexi), lexi_to_num(b_bwd_lexi))
        num_min = min(lexi_to_num(a_bwd_lexi), lexi_to_num(b_bwd_lexi))
        b_words = Vector{Int}[]
        # iterate over all binary strings of length bwd_len that are between min_bwd_lexi and max_bwd_lexi in lexicographical order
        for i in num_min:num_max
            push!(b_words, num_to_lexi(i, bwd_len) |> _lexi_to_knead)
        end
        # amalgamate pairs of words that differ only in the first symbol, since they represent the same pruning condition
        for _ in 1:bwd_len
            tmp_words = Vector{Int}[]
            for w in b_words
                w_alt = copy(w)
                w_alt[end] = 1 - w_alt[end]
                if w_alt in b_words
                    push!(tmp_words, w_alt[1:end-1])
                    # remove w and w_alt from b_words to avoid duplication
                    deleteat!(b_words, findfirst(==(w), b_words))
                    deleteat!(b_words, findfirst(==(w_alt), b_words))
                else
                    push!(tmp_words, w)
                end
            end
            b_words = tmp_words
        end

        # compute forward part
        fwd_len = max(length(b.forward), length(c.forward))
        b_fwd_lexi = _knead_to_lexi(b[0:fwd_len-1])
        c_fwd_lexi = _knead_to_lexi(c[0:fwd_len-1])
        num_max = max(lexi_to_num(b_fwd_lexi), lexi_to_num(c_fwd_lexi))
        num_min = min(lexi_to_num(b_fwd_lexi), lexi_to_num(c_fwd_lexi))
        f_words = Vector{Int}[]
        for i in num_min:num_max
            push!(f_words, num_to_lexi(i, fwd_len) |> _lexi_to_knead)
        end
        # amalgamate pairs of words that differ only in the first symbol, since they represent the same pruning condition
        for _ in 1:fwd_len
            tmp_words = Vector{Int}[]
            for w in f_words
                w_alt = copy(w)
                w_alt[end] = 1 - w_alt[end]
                if w_alt in f_words
                    push!(tmp_words, w_alt[1:end-1])
                    deleteat!(f_words, findfirst(==(w), f_words))
                    deleteat!(f_words, findfirst(==(w_alt), f_words))
                else
                    push!(tmp_words, w)
                end
            end
            f_words = tmp_words
        end
        # combine backward and forward parts to get forbidden words
        for bwd in reverse.(b_words)
            for fwd in f_words
                push!(words, string(join(bwd, ""), join(fwd, "")))
            end
        end
    end
    
    return words
end

"""
    markerize(words::Vector{String}) -> Vector{String}

Given a list of forbidden words, compute the corresponding markers.

# Example
```jldoctest pruning
julia> markerize(["0010100", "0011100", "0010101"])
2-element Vector{String}:
 "001010X"
 "0011100"
```
"""
function markerize(words::Vector{String})
    markers = Vector{String}()
    while !isempty(words)
        w = pop!(words)
        found_pair = false
        for i in 1:length(w)
            target = collect(w)
            target[i] = w[i] == '0' ? '1' : '0'
            if String(target) in words
                push!(markers, w[1:i-1] * "X" * w[i+1:end])
                deleteat!(words, findfirst(==(String(target)), words))
                found_pair = true
                break
            end
        end
        found_pair || push!(markers, w)
    end
    return markers
end
