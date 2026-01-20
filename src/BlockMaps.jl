"`BlockMaps.jl` provides a library of block maps (sliding block code) on subshift of finite type (SFT) aiming for easy-to-use and fast implementation."
module BlockMaps
# using some packages
import Base: length, inv, Matrix, oneunit
using LinearAlgebra: pinv, tr, I
using Permutations

export VertexWalks, EdgeWalks, VertexCycles, EdgeCycles
export BlockMap
export marker_endomorphism, shift
export shortenable, shorten, extend
export ElementaryStrongShiftEquivalence, ESSE
export StrongShiftEquivalence, SSE
export reverse, reverse!, reduce, reduce!, reducefirst, reducefirst!, reducelast, reducelast!
export dimension_representation
export orbitsign_number, os, gyration_number, gy, sign_gyration_compatibility_condition, sgcc

include("matrix2str.jl")

# TODO: ファイルを適切に分割
# TODO: テストを書く
# TODO: ドキュメントを書く
# TODO: パッケージ化

"abstract iterator type to walk on graph"
abstract type GraphWalks end
abstract type CycleWalks end

"iterator type to walk on vertices in graph"
struct VertexWalks <: GraphWalks
    A :: Matrix{Int} # Adjacency matrix
    n :: Int # length of the walk
    start :: Vector{Int} # the first vertex sequence
    function VertexWalks(A::Matrix{Int}, n::Int, start::Vector{Int})
        # validation check
        if size(A)[1] != size(A)[2]
            throw(ArgumentError("adjacency matrix must be square"))
        end
        if n < 1
            throw(ArgumentError("the length of walks must be positive"))
        end
        if length(start) == 0 && sum(A^n) == 0
            return new(A, n, start)
        elseif length(start) != n
            throw(ArgumentError("length of start sequence is not equal to $n"))
        end
        if all(1 .<= start .<= size(A)[1]) |> !
            throw(ArgumentError("each symbol must be between 1 and size of graph"))
        end
        # check start seq is valid
        for i in 1:length(start)-1
            if A[start[i], start[i+1]] == 0
                throw(ArgumentError("start sequence must be valid in graph"))
            end
        end
        return new(A, n, start)
    end
end

function VertexWalks(A::Matrix{Int}, n::Int, start::Tuple{Vararg{Int}})
    return VertexWalks(A, n, [start...])
end

# Constructor of VertexWalks without start sequence
function VertexWalks(A::Matrix{Int}, n::Int)
    walk = Vector{Int}(undef, n)
    N = size(A)[1]
    for s in 1:n
        walk[1] = s
        found_walk = true
        for i in 1:n-1
            try
                walk[i+1] = findfirst(x -> x > 0, A[walk[i], :])
            catch
                found_walk = false
                break
            end
        end
        if found_walk return VertexWalks(A, n, walk) end
    end
    return VertexWalks(A, n, (Vector{Int}(undef, 0)))
end

function Base.iterate(iter::VertexWalks)
    if sum(iter.A ^ iter.n) == 0 return nothing end
    return tuple(iter.start...), iter.start
end

function Base.iterate(iter::VertexWalks, state::Vector{Int})
    walk = state
    N = size(iter.A)[1]
    if iter.n == 1 && walk[1] >= N
        return nothing
    elseif iter.n == 1
        return tuple(walk[1] + 1), walk .+ 1
    end

    # try increment the walk
    for d in reverse(2:iter.n)
        # try increment d th vertex
        prev_dth_vertex = walk[d]
        tmp_d = findfirst(x -> x > 0, iter.A[walk[d-1], walk[d]+1:end])
        if !isnothing(tmp_d)
            walk[d] = tmp_d + prev_dth_vertex
            # set the next vertex to the smallest one
            for i in d+1:iter.n
                walk[i] = findfirst(x -> x > 0, iter.A[walk[i-1], :])
            end
            return tuple(walk...), walk
        end
    end
    for _ in walk[1]:N-1
        walk[1] += 1
        found_walk = true
        for d in 2:iter.n
            tmp_d = findfirst(x -> x > 0, iter.A[walk[d-1], :])
            if isnothing(tmp_d)
                found_walk = false
                break
            end
            walk[d] = tmp_d
        end
        if found_walk return tuple(walk...), walk end
    end

    return nothing # end of iteration
end

"Compute indices of start and end vertex in given edge index"
function _get_vertex_index(A::Matrix{Int}, e::Int)
    if !(1 <= e <= sum(A)) throw(ArgumentError("edge is not in graph")) end
    s = 0
    local ind
    for (i, a) in enumerate(A)
        if e <= s + a
            ind = i
            break
        end
        s += a
    end
    v_s = ((ind - 1) % size(A)[1]) + 1
    v_t = ((ind - 1) ÷ size(A)[1]) + 1
    return v_s, v_t
end

"Iterator type to walk on edges in graph"
struct EdgeWalks <: GraphWalks
    A :: Matrix{Int} # adjacency matrix
    n :: Int # length of each walk
    start :: Vector{Int} #  the first edge sequence

    function EdgeWalks(A::Matrix{Int}, n::Int, start::Vector{Int})
        # validation check
        if size(A)[1] != size(A)[2]
            throw(ArgumentError("adjacency matrix must be square"))
        end
        if length(start) <= 0 && sum(A^n) <= 0
            new(A, n, start)
        end
        if n < 1
            throw(ArgumentError("the length of the walk must be positive"))
        end
        if length(start) != n
            throw(ArgumentError("length of start sequence is not equal to $n"))
        end
        if all(1 .<= start .<= sum(A)) |> !
            throw(ArgumentError("each symbol must be between 1 and sum of adjacency matrix"))
        end
        # check start is valid
        _, end_vertex = _get_vertex_index(A, start[1])
        for i in 1:n-1
            start_vertex, tmp = _get_vertex_index(A, start[i+1])
            if end_vertex != start_vertex
                throw(ArgumentError("start sequence is not valid in graph"))
            end
            end_vertex = tmp
        end
        new(A, n, start)
    end
end

function _edge_index(A::Matrix)
    to_return::Vector{Union{UnitRange{Int}, Nothing}} = [1+i:j for (i, j) in zip(cumsum([0, A...][1:end-1]), cumsum([A...]))]
    to_return[[A...] .<= 0] .= nothing
    return reshape(to_return, size(A))
end

function EdgeWalks(A::Matrix{Int}, n::Int, start::Tuple{Vararg{Int}})
    return EdgeWalks(A, n, [start...])
end

# Constructor of EdgeWalks without start sequence
function EdgeWalks(A::Matrix{Int}, n::Int)
    walk = Vector{Int}(undef, n)
    S = sum(A)
    edge_index = _edge_index(A)
    for s in 1:S
        walk[1] = s
        found_walk = true
        for i in 1:n-1
            try
                walk[i+1] = findfirst(x -> _get_vertex_index(A, walk[i])[2] == _get_vertex_index(A, x)[1], 1:S)
            catch
                found_walk = false
                break
            end
        end
        if found_walk return EdgeWalks(A, n, walk) end
    end
    return EdgeWalks(A, n, Vector{Int}(undef, 0))
end

function Base.iterate(iter::EdgeWalks)
    if sum(iter.A ^ iter.n) == 0 return nothing end
    edge_index = _edge_index(iter.A)
    reshape(edge_index, size(iter.A))
    return tuple(iter.start...), (iter.start, edge_index)
end

function Base.iterate(iter::EdgeWalks, state)
    walk::Vector{Int} = state[1]
    edge_index::Matrix = state[2]
    S = sum(iter.A)
    if iter.n == 1 && walk[1] >= S
        return nothing
    elseif iter.n == 1
        return tuple(walk[1] + 1), (walk .+ 1, edge_index)
    end
    for d in reverse(2:iter.n)
        prev_dth_edge = walk[d]
        v_s, v_e = _get_vertex_index(iter.A, prev_dth_edge)
        if prev_dth_edge+1 in edge_index[v_s, v_e]
            walk[d] = prev_dth_edge + 1
            for i in d+1:iter.n
                _, v_e = _get_vertex_index(iter.A, walk[i-1])
                walk[i] = findfirst(x -> v_e == _get_vertex_index(iter.A, x)[1], 1:S)
            end
            return tuple(walk...), (walk, edge_index)
        else
            if v_e + 1 > size(iter.A)[1]
                continue
            end
            tmp_d = findfirst(!isnothing, edge_index[v_s, v_e+1:end])
            if isnothing(tmp_d) continue end
            walk[d] = first(edge_index[v_s, v_e+tmp_d])
            for i in d+1:iter.n
                _, v_e = _get_vertex_index(iter.A, walk[i-1])
                walk[i] = findfirst(x -> v_e == _get_vertex_index(iter.A, x)[1], 1:S)
            end
            return tuple(walk...), (walk, edge_index)
        end
    end
    for _ in walk[1]:S-1
        walk[1] += 1
        _, v_e = _get_vertex_index(iter.A, walk[1])
        found_walk = true
        for d in 2:iter.n
            try
                walk[d] = findfirst(x -> v_e == _get_vertex_index(iter.A, x)[1], 1:S)
            catch
                found_walk = false
                break
            end
            _, v_e = _get_vertex_index(iter.A, walk[d])
        end
        if found_walk return tuple(walk...), (walk, edge_index) end
    end
    # end of iterate
    return nothing
end

Base.length(iter::GraphWalks) = sum(iter.A ^ iter.n)
Base.eltype(iter::GraphWalks) = NTuple{iter.n, Int}


function _is_period(s::Tuple{Vararg{Int}}, n::Int)
    (length(s) % n != 0) && return false
    for i in 1:length(s)÷n
        if s[1:n] != s[(i-1)*n+1:i*n]
            return false
        end
    end

    return true
end

function _least_periodof(s::Tuple{Vararg{Int}})
    n = length(s)
    for i in 1:n
        if _is_period(s, i)
            return i
        end
    end
    # not periodic case
    return 0
end

struct VertexCycles <: CycleWalks
    A :: Matrix{Int} # adjacency matrix
    n :: Int # length of each walk
    start :: Vector{Int} #  the first vertex sequence
    least :: Bool # whether the period is the least or not

    function VertexCycles(A::Matrix{Int}, n::Int, start::Vector{Int}, least::Bool=true)
        # validation check
        if size(A)[1] != size(A)[2]
            throw(ArgumentError("Adjacency matrix must be square"))
        end

        _l = tr(A^n)
        for i in 1:(n-1)
            _l -= n % i == 0 ? tr(A^i) : 0
        end
        if length(start) <= 0 && _l <= 0
            new(A, n, start)
        end

        if n < 1
            throw(ArgumentError("The length of cycles must be positive"))
        end
        if length(start) != n
            throw(ArgumentError("length of start walk is not equal to $n"))
        end
        if all(1 .<= start .<= size(A)[1]) |> !
            throw(ArgumentError("each symbol must be between 1 and size of adjacency matrix"))
        end
        # check start seq is valid
        for i in 1:length(start)-1
            if A[start[i], start[i+1]] == 0
                throw(ArgumentError("start sequence must be valid in graph"))
            end
        end
        if A[start[end], start[1]] == 0
            throw(ArgumentError("start sequence must be periodic in graph"))
        end
        if least && _least_periodof(start) != n
            throw(ArgumentError("least period of start sequence must be the $n"))
        end
        new(A, n, start, least)
    end
end

# Constructor of VertexCycles without start sequence
function VertexCycles(A::Matrix{Int}, n::Int, least::Bool=true)
    for w in VertexWalks(A, n)
        if A[w[end], w[1]] != 0
            if (least && _least_periodof(w) == n) || !least
                return VertexCycles(A, n, [w...], least)
            end
        end
    end
    return VertexCycles(A, n, Vector{Int}(undef, 0), least)
end

function Base.iterate(iter::VertexCycles)
    if length(iter) == 0 return nothing end
    return tuple(iter.start...), iter.start
end

function Base.iterate(iter::VertexCycles, state::Vector{Int})
    for (i, w) in enumerate(VertexWalks(iter.A, iter.n, state))
        if i > 1 && iter.A[w[end], w[1]] != 0
            return w, [w...,]
        end
    end

    return nothing
end

struct EdgeCycles <: CycleWalks
    A :: Matrix{Int} # adjacency matrix
    n :: Int # length of each walk
    start :: Vector{Int} #  the first edge sequence
    least :: Bool # whether the period is the least or not

    function EdgeCycles(A::Matrix{Int}, n::Int, start::Vector{Int}, least=true)
        # validation check
        if size(A)[1] != size(A)[2]
            throw(ArgumentError("Adjacency matrix must be square"))
        end

        _l = tr(A^n)
        for i in 1:(n-1)
            _l -= n % i == 0 ? tr(A^i) : 0
        end
        if length(start) <= 0 && _l <= 0
            new(A, n, start)
        end
        
        if n < 1
            throw(ArgumentError("The length of cycles must be positive"))
        end
        if length(start) != n
            throw(ArgumentError("length of start walk is not equal to $n"))
        end
        if all(1 .<= start .<= sum(A)) |> !
            throw(ArgumentError("each symbol must be between 1 and sum of adjacency matrix"))
        end
        # check start seq is valid
        _, end_vertex = _get_vertex_index(A, start[1])
        for i in 1:n-1
            start_vertex, tmp = _get_vertex_index(A, start[i+1])
            if end_vertex != start_vertex
                throw(ArgumentError("start sequence is not valid in graph"))
            end
            end_vertex = tmp
        end
        if end_vertex != _get_vertex_index(A, start[1])[1]
            throw(ArgumentError("start sequence must be periodic in graph"))
        end
        if least && _least_periodof(tuple(start...)) != n
            throw(ArgumentError("least period of start sequence must be the $n"))
        end
        new(A, n, start, least)
    end
end

# Constructor of EdgeCycles without start sequence
function EdgeCycles(A::Matrix{Int}, n::Int, least=true)
    for w in EdgeWalks(A, n)
        if _get_vertex_index(A, w[end])[2] == _get_vertex_index(A, w[1])[1]
            if (least && _least_periodof(w) == n) || !least
                return EdgeCycles(A, n, [w...], least)
            end
        end
    end
    return EdgeCycles(A, n, Vector{Int}(undef, 0), least)
end

function Base.iterate(iter::EdgeCycles)
    if length(iter) == 0 return nothing end
    return tuple(iter.start...), iter.start
end

function Base.iterate(iter::EdgeCycles, state)
    for (i, w) in enumerate(EdgeWalks(iter.A, iter.n, state))
        if i > 1 && _get_vertex_index(iter.A, w[end])[2] == _get_vertex_index(iter.A, w[1])[1]
            if iter.least && _least_periodof(w) != iter.n
                continue
            end
            return w, [w...]
        end
    end

    return nothing
end

function Base.length(iter::CycleWalks)
    if iter.least
        l = tr(iter.A ^ iter.n)
        for i in 1:(iter.n-1)
            if l % i == 0
                l -= tr(iter.A ^ i)
            end
        end
        return l
    else
        return tr(iter.A ^ iter.n)
    end
end
Base.eltype(iter::CycleWalks) = NTuple{iter.n, Int}

"""
    BlockMap{N}(start::Matrix{Int}, target::Matrix{Int}, blocks::Dict{NTuple{N, Int}, Int}, memory::Int)

Construct a N-block map from the given start matrix to the target matrix with the given blocks.
If matrix contain only 0 and 1, we treat SFT as vertex shift. If not, we treat SFT as edge shift.
For edge shift, the index of the edge is ordered from top to bottom and left to right in the matrix (column major).
"""
struct BlockMap{N}
    start :: Matrix{Int} # definition domain of the map
    target :: Matrix{Int} # codomain of the map
    blocks :: Dict{NTuple{N, Int}, Int}
    memory :: Int # how far refer the past symbol to determine the current symbol
    function BlockMap(start::Matrix{Int}, target::Matrix{Int}, blocks::Dict{NTuple{N, Int}, Int}, memory::Int = 0) where N
        _check_block_map(start, target, blocks)
        return new{N}(start, target, blocks, memory)
    end
end

function _check_vertex_map_input(start, target, blocks)
    A = start
    n = size(A)[1]
    l = length(collect(keys(blocks))[1])

    # check that keys has just valid words of A and output sequence is valid
    for w in Iterators.product(fill(1:n, l)...)
        # check keys has w if w is a valid word of A
        is_valid = true
        for i in 1:length(w)-1
            if A[w[i], w[i+1]] == 0
                is_valid = false
                break
            end
        end
        if is_valid
            if !haskey(blocks, w)
                throw(ArgumentError("block map must have all valid words of the start shift space"))
            end
        else
            if haskey(blocks, w)
                throw(ArgumentError("block map must not have invalid words of the start shift space"))
            end
        end
    end

    return true
end

function _check_edge_map_input(start, target, blocks)
    A = start
    n = size(A)[1]
    l = length(collect(keys(blocks))[1])

    _cumsumA = cumsum([0, A...])
    ind2edge_range = [_cumsumA[i]+1:s for (i, s) in enumerate(cumsum([A...]))]

    _checked = Vector{Tuple{Vararg{Int}}}()
    for w_v in Iterators.product(fill(1:n, l+1)...)
        w_e = [ind2edge_range[v + n * (w_v[i+1] - 1)] for (i, v) in enumerate(w_v[1:end-1])]
        for e_seq in Iterators.product(w_e...)
            if !haskey(blocks, e_seq)
                throw(ArgumentError("block map must have all valid words of the start shift space"))
            end
            push!(_checked, e_seq)
        end
    end
    if length(_checked) != length(blocks)
        throw(ArgumentError("block map must not have invalid words of the start shift space"))
    end
end

function _check_vertex_map_output(start, target, blocks)
    A = start
    B = target
    m = size(B)[1]
    n = size(A)[1]
    l = length(collect(keys(blocks))[1])

    for (x1, y1) in blocks, (x2, y2) in blocks
        if !((1 <= y1 <= m) && (1 <= y2 <= m))
            throw(ArgumentError("block map must output numbers from 1 to the size of the target matrix"))
        end
        if (l > 1) && (x1[2:end] == x2[1:end-1])
            if B[y1, y2] == 0
                throw(ArgumentError("block map must not output invalid sequence in the target shift space"))
            end
        elseif l == 1
            if _is_on_vertex_shift(start, target, blocks)[1]
                if A[x1[1], x2[1]] == 1 && B[y1, y2] == 0
                    throw(ArgumentError("block map must not output invalid sequence in the target shift space"))
                end
            else
                _cumsumA = cumsum([0, A...])
                ind2edge_range = [_cumsumA[i]+1:s for (i, s) in enumerate(cumsum([A...]))]
                u = findfirst(r -> x1[1] in r, ind2edge_range)
                u = ((u - 1) ÷ n) + 1
                v = findfirst(r -> x2[1] in r, ind2edge_range)
                v = ((v - 1) % n) + 1
                if u == v && B[y1, y2] == 0
                    throw(ArgumentError("block map must not output invalid sequence in the target shift space"))
                end
            end
        end
    end

    return true
end

function _check_edge_map_output(start, target, blocks)
    A = start
    B = target
    m = size(B)[1]
    n = size(A)[1]
    l = length(collect(keys(blocks))[1])

    _cumsumB = cumsum([0, B...])
    ind2edge_range = [_cumsumB[i]+1:s for (i, s) in enumerate(cumsum([B...]))]

    for (x1, y1) in blocks, (x2, y2) in blocks
        if !(1 <= y1 <= sum(B)) || !(1 <= y2 <= sum(B))
            throw(ArgumentError("block map must output a number between 1 to the sum of target matrix"))
        end
        u = findfirst(r -> y1 in r, ind2edge_range)
        u = ((u - 1) ÷ m) + 1 # end vertex of y1
        v = findfirst(r -> y2 in r, ind2edge_range)
        v = ((v - 1) % m) + 1 # start vertex of y2
        if ((l > 1) && (x1[2:end] == x2[1:end-1])) && (u != v)
            throw(ArgumentError("block map must not output invalid sequence in target shift space"))
        elseif (l == 1) && (u != v)
            if _is_on_vertex_shift(start, target, blocks)[1]
                if A[x1[1], x[2]] == 1
                    throw(ArgumentError("block map must not output invalid sequence in target shift space"))
                end
            else
                _cumsumA = cumsum([0, A...])
                ind2edge_rangeA = [_cumsumA[i]+1:s for (i, s) in enumerate(cumsum([A...]))]
                u_x = findfirst(r -> x1[1] in r, ind2edge_rangeA)
                u_x = ((u_x - 1) ÷ n) + 1
                v_x = findfirst(r -> x2[1] in r, ind2edge_rangeA)
                v_x = ((v_x - 1) % n) + 1
                if u_x == v_x
                    throw(ArgumentError("block map must not output invalid sequence in target shift space"))
                end
            end
        end
    end

    return true
end

"throw errors if the block_map is not well-defined"
function _check_block_map(start::Matrix{Int}, target::Matrix{Int}, blocks::Dict{NTuple{N, Int}, Int}) where N
    # check matrix start and target is square
    if size(start)[1] != size(start)[2]
        throw(ArgumentError("start matrix must be a square"))
    end
    if size(target)[1] != size(target)[2]
        throw(ArgumentError("target matrix must be a square"))
    end

    # check the all blocks have the same length
    if length(Set(length.(keys(blocks)))) > 1
        throw(ArgumentError("all blocks must have the same length"))
    end

    is_v = _is_on_vertex_shift(start, target, blocks)

    if is_v[1]
        _check_vertex_map_input(start, target, blocks)
    else
        _check_edge_map_input(start, target, blocks)
    end

    if is_v[2]
        _check_vertex_map_output(start, target, blocks)
    else
        _check_edge_map_output(start, target, blocks)
    end

    return true
end

function Base.show(io::IO, f::BlockMap)
    print(io, "BlockMap(start=$(f.start), target=$(f.target), (")

    B = collect(f.blocks) |> sort
    for xy in B[1:end-1]
        print(io, "$xy, ")
    end
    print(io, B[end])

    print(io, "), memory=$(f.memory))")
end

function Base.show(io::IO, ::MIME"text/plain", f::BlockMap)
    println(io, "$(length(f))-BlockMap: SFT($(f.start)) → SFT($(f.target)) with memory=$(f.memory):")

    d_row = displaysize(stdout)[1] - 5 # -5 is for input prompt, header, ⋮, empty line and new input line

    B = collect(f.blocks) |> sort
    if length(B) < d_row
        for xy in B[1:end-1]
            println(io, xy)
        end
        print(io, B[end])
    else
        for xy in B[1:(d_row+1)÷2]
            println(io, xy)
        end
        l = length(string(B[1]))
        println(io, repeat(" ", l÷2) * "⋮")
        for xy in B[end-((d_row+1)÷2)+1:end-1]
            println(io, xy)
        end
        print(io, B[end])
    end
end

"""
    _is_on_vertex_shift(start, target, blocks)

Judge whether the given block map is a vertex shift or not

Return: `(start is vertex shift, target is vertex shift)`
"""
function _is_on_vertex_shift(start, target, blocks)
    for_start = true
    for_target = true
    if all(0 .<= start .<= 1)
    else
        for_start = false
    end
    if all(0 .<= target .<= 1)
    else
        for_target = false
    end

    n = size(start)[1]
    m = size(target)[1]

    for x in collect(keys(blocks))
        for c in x
            if !(1 <= c <= n)
                for_start = false
            end
        end
    end
    for y in collect(values(blocks))
        if !(1 <= y <= m)
            for_target = false
        end
    end

    return (for_start, for_target)
end

_is_on_vertex_shift(block_map::BlockMap) = _is_on_vertex_shift(block_map.start, block_map.target, block_map.blocks)

"""
    BlockMap(A::Matrix{Int}, B::Matrix{Int}, blocks::Dict{String, Int}, memory::Int)

Construct a block map from symbolic strings as input blocks and integers as output symbols. 
The string has numbers as symbol and each symbol separated by space.

# Examples
```jldoctest
julia> BlockMap(
        [2;;], 
        [2;;], 
        Dict("1 1" => 1, "1 2" => 2, "2 1" => 1, "2 2" => 2), 
        0
    )
2-BlockMap: SFT([2;;]) → SFT([2;;]) with memory=0:
(1, 1) => 1
(1, 2) => 2
(2, 1) => 1
(2, 2) => 2
"""
function BlockMap(A::Matrix{Int}, B::Matrix{Int}, blocks::Dict{String, Int}, memory::Int)
    if split.(collect(keys(blocks)), " ") |> length |> unique |> length > 1
        throw(ArgumentError("all blocks must have the same length"))
    end

    N = split(first(keys(blocks)), " ") |> length

    proper_blocks = Dict{NTuple{N, Int}, Int}()
    for (k, v) in blocks
        k = parse.(Int, split(k, " "))
        k = (k...,)
        proper_blocks[k] = v
    end
    ϕ = BlockMap(A, B, proper_blocks, memory)
    return ϕ
end

"alias of defining blocks with Char"
function BlockMap(A::Matrix{Int}, B::Matrix{Int}, blocks::Dict{String, Char}, memory::Int)
    if split.(collect(keys(blocks)), " ") |> length |> unique |> length > 1
        throw(ArgumentError("all blocks must have the same length"))
    end

    N = split(first(keys(blocks)), " ") |> length

    proper_blocks = Dict{NTuple{N, Int}, Int}()
    for (k, v) in blocks
        k = parse.(Int, split(k, " "))
        k = (k...,)
        proper_blocks[k] = parse(Int, v)
    end
    ϕ = BlockMap(A, B, proper_blocks, memory)
    return ϕ
end

# Constructor of endo-block map
BlockMap(A::Matrix{Int}, blocks, memory::Int) = BlockMap(A, A, blocks, memory)

Base.length(block_map::BlockMap{N}) where N = N

"""
    (block_map::BlockMap)(seq::Int...) -> ::Tuple{Vararg{Int}}

Return the block map output of given symbolic sequence.

# Examples
```jldoctest
```
"""
function (block_map::BlockMap)(in_seq::Tuple{Vararg{Int}})::Tuple{Vararg{Int}}
    l = length(block_map)
    out_len = length(in_seq) - l + 1
    out_seq = Vector{Int}(undef, out_len)
    for i in 1:out_len
        block = in_seq[i:i+l-1]
        if haskey(block_map.blocks, block)
            out_seq[i] = block_map.blocks[block]
        else
            throw(ArgumentError("The input sequence is not valid"))
        end
    end

    return (out_seq...,)
end

function (block_map::BlockMap)(in_seq::Int...)
    return block_map(in_seq)
end

function (block_map::BlockMap)(in_seq::Vector{Int})
    return block_map(in_seq...)
end

function (block_map::BlockMap)(in_seq::String)
    return block_map(parse.(Int, split(in_seq, " "))...)
end

function (block_map::BlockMap)(in_seq::Dict{Int, Int})
    in_seq = last.(sort(collect(in_seq), by=first))
    out_seq = block_map(collect(values(in_seq))...)
    first_ind = minimum(keys(in_seq))
    to_return = Dict{Int, Int}()
    for (i, y) in enumerate(out_seq)
        to_return[first_ind + i - 1 + block_map.memory] = y
    end
    return to_return
end


function _shift_orbits(iter::CycleWalks)
    orbits = Vector{Vector{NTuple{iter.n, Int}}}()
    for w in iter
        nw = [w...]
        find_orb = false
        for (i, o) in enumerate(orbits), _ in 1:iter.n
            circshift!(nw, 1)
            if tuple(nw...) in o
                push!(orbits[i], w)
                find_orb = true
                break
            end
        end
        !find_orb && push!(orbits, [w])
    end

    return orbits
end

function vertex_orbits(A::Matrix{Int}, n::Int, least::Bool=true)
    return _shift_orbits(VertexCycles(A, n, least))
end

function edge_orbits(A::Matrix{Int}, n::Int, least::Bool=true)
    return _shift_orbits(EdgeCycles(A, n, least))
end

"""
    shift_orbits(f::BlockMap, n::Int, least::Bool=true)

Compute the shift n-orbits in the SFT of the given block map domain/codomain.
"""
function shift_orbits(f::BlockMap, n::Int, least::Bool=true)

    is_start_vertex, is_target_vertex = _is_on_vertex_shift(f)

    local start_Os, target_Os
    if is_start_vertex
        start_Os = vertex_orbits(f.start, n, least)
    else
        start_Os = edge_orbits(f.start, n, least)
    end
    if is_target_vertex
        target_Os = vertex_orbits(f.target, n, least)
    else
        target_Os = edge_orbits(f.target, n, least)
    end

    return start_Os, target_Os
end

"""
    ElementaryStrongShiftEquivalence

Elementary strong shift equivalence is a pair of matrices.
"""
mutable struct ElementaryStrongShiftEquivalence{T}
    first :: Matrix{T}
    second :: Matrix{T}

    function ElementaryStrongShiftEquivalence{T}(esse::Pair{Matrix{T}, Matrix{T}}) where T
        R, S = esse
        # validation
        if size(R)[1] != size(S)[2]
            throw(ArgumentError("The matrices must be able to multiple of each other"))
        elseif size(S)[1] != size(R)[2]
            throw(ArgumentError("The matrices must be able to multiple of each other"))
        else
            return new{T}(R, S)
        end
    end
end
ElementaryStrongShiftEquivalence(esse::Pair{Matrix{T}, Matrix{T}}) where T = ElementaryStrongShiftEquivalence{T}(esse)
ElementaryStrongShiftEquivalence(R::Matrix{T}, S::Matrix{T}) where T = ElementaryStrongShiftEquivalence{T}(R => S)
"Alias of `ElementaryStrongShiftEquivalence`"
const ESSE{T} = ElementaryStrongShiftEquivalence{T}
Base.iterate(esse::ESSE, i=1) = i > 2 ? nothing : (getfield(esse, i), i + 1)
Base.indexed_iterate(esse::ESSE, i::Int, state=1) = (getfield(esse, i), i + 1)
Base.length(esse::ESSE) = 2
Base.getindex(esse::ESSE, i::Int) = (esse.first => esse.second)[i]
Base.first(esse::ESSE) = esse.first
Base.last(esse::ESSE) = esse.second
Base.firstindex(esse::ESSE) = 1
Base.lastindex(esse::ESSE) = 2
Base.convert(::Type{ESSE{T}}, esse::Pair{Matrix{T}, Matrix{T}}) where T = ESSE(esse)
function Base.reverse!(esse::ESSE)
    tmp = esse.first
    esse.first = esse.second
    esse.second = tmp
    return esse
end
Base.reverse(esse::ESSE) = ESSE(esse.second => esse.first)
Base.convert(::Type{Pair{Matrix{T}, Matrix{T}}}, esse::ESSE{T}) where T = (esse.first => esse.second)

function Base.show(io::IO, esse::ESSE)
    print(io, "ElementaryStrongShiftEquivalence(")
    print(io, esse.first, " => ", esse.second)
    print(io, ")")
end

"""
    show(io::IO, ::MIME"text/plain", esse::ESSE)

Print the given ElementaryStrongShiftEquivalence object.

# Examples
```jldoctest
julia> R = [1 1; 0 1]
2×2 Matrix{Int64}:
 1  1
 0  1

julia> S = [1 0; 1 1]
2×2 Matrix{Int64}:
 1  0
 1  1

julia> esse = ESSE(R, S)
ElementaryStrongShiftEquivalence:
      ⎡1 1⎤ ⎡1 0⎤
⎡2 1⎤ ⎣0 1⎦ ⎣1 1⎦  ⎡1 1⎤
⎣1 1⎦─────────────▶︎⎣1 2⎦
"""
function Base.show(io::IO, ::MIME"text/plain", esse::ESSE)
    d_h, d_w = displaysize(stdout)

    println(io, "ElementaryStrongShiftEquivalence:")
    R, S = esse
    RS, SR = R * S, S * R

    R_str = split(matrix2str(R), "\n")
    S_str = split(matrix2str(S), "\n")
    RS_str = split(matrix2str(RS), "\n")
    SR_str = split(matrix2str(SR), "\n")

    h = max(length(R_str)+1, length(S_str)+1)
    w = length(R_str[1]) + length(S_str[1]) + length(RS_str[1]) + length(SR_str[1]) + 4

    show_with_braille = false

    if h + 4 > d_h
        show_with_braille = true
    elseif w > d_w
        show_with_braille = true
    end

    RS_str = show_with_braille ? split(matrix2braille(RS), "\n") : RS_str
    R_str = show_with_braille ? split(matrix2braille(R), "\n") : R_str
    S_str = show_with_braille ? split(matrix2braille(S), "\n") : S_str
    SR_str = show_with_braille ? split(matrix2braille(SR), "\n") : SR_str

    h = max(length(R_str)+1, length(S_str)+1)
    w = length(R_str[1]) + length(S_str[1]) + length(RS_str[1]) + length(SR_str[1]) + 4

    s = 1.0
    if h + 4 > d_h || w > d_h
        s = min(d_h / (h+3), d_w / w)
    end

    RS_str = show_with_braille ? split(matrix2braille(RS, floor(Int, s * length(RS_str[1])), floor(Int, s * length(RS_str))), "\n") : RS_str
    R_str = show_with_braille ? split(matrix2braille(R, floor(Int, s * length(R_str[1])), floor(Int, s * length(R_str))), "\n") : R_str
    S_str = show_with_braille ? split(matrix2braille(S, floor(Int, s * length(S_str[1])), floor(Int, s * length(S_str))), "\n") : S_str
    SR_str = show_with_braille ? split(matrix2braille(SR, floor(Int, s * length(SR_str[1])), floor(Int,s * length(SR_str))), "\n") : SR_str

    max_h = max(length(R_str)+1, length(S_str)+1)
    for i in reverse(2:max_h)
        if i <= length(RS_str)
            j = length(RS_str) - i + 1
            print(io, RS_str[j])
        else
            print(io, " "^(length(RS_str[1])))
        end
        print(io, " ")
        if i <= length(R_str) + 1
            j = length(R_str) - i + 2
            print(io, R_str[j])
        else
            print(io, " "^(length(R_str[1])))
        end
        print(io, " ")
        if i <= length(S_str) + 1
            j = length(S_str) - i + 2
            print(io, S_str[j])
        else
            print(io, " "^(length(S_str[1])))
        end
        print(io, "  ")
        if i <= length(SR_str)
            j = length(SR_str) - i + 1
            println(io, SR_str[j])
        else
            println(io, " "^(length(SR_str[1])))
        end
    end
    print(io, RS_str[end])
    print(io, "─"^(length(R_str[1]) + length(S_str[1]) + 3))
    print(io, "▶︎", SR_str[end])
end

"""
    StrongShiftEquivalence

Strong shift equivalence is a sequence of elementary strong shift equivalences.
"""
mutable struct StrongShiftEquivalence{T}
    sequence::Vector{ESSE{T}}

    function StrongShiftEquivalence{T}(sequence::ESSE{T}...) where T
        # validation
        X = sequence[1][2] * sequence[1][1]
        for (i, (R, S)) in enumerate(sequence[2:end])
            if X != R * S
                throw(ArgumentError("$(i)th and $(i+1)th matrix don't satisfy the strong shift equations"))
            end
            X = S * R
        end
        return new{T}([sequence...])
    end
end
"Alias of StrongShiftEquivalence"
const SSE{T} = StrongShiftEquivalence{T} # alias
StrongShiftEquivalence(sequence::ESSE{T}...) where T = StrongShiftEquivalence{T}(sequence...)
StrongShiftEquivalence(sequence::Vector{ESSE{T}}) where T = StrongShiftEquivalence{T}(sequence...)
Base.iterate(sse::SSE, args...) = iterate(sse.sequence, args...)
Base.length(sse::SSE) = length(sse.sequence)
Base.getindex(sse::SSE, i::Int) = getindex(sse.sequence, i)
Base.first(sse::SSE) = first(sse.sequence)
Base.last(sse::SSE) = last(sse.sequence)
Base.firstindex(sse::SSE) = firstindex(sse.sequence)
Base.lastindex(sse::SSE) = lastindex(sse.sequence)
Base.convert(::Type{SSE{T}}, sse::Vector{ESSE{T}}) where T = SSE(sse)
Base.convert(::Type{Vector{ESSE{T}}}, sse::SSE{T}) where T = sse.sequence
Base.pop!(sse::SSE) = pop!(sse.sequence)
Base.popfirst!(sse::SSE) = popfirst!(sse.sequence)
Base.popat!(sse::SSE, i::Int) = popat!(sse.sequence, i)
Base.push!(sse::SSE, esse::ESSE) = push!(sse.sequence, esse)
Base.pushfirst!(sse::SSE, esse::ESSE) = pushfirst!(sse.sequence, esse)
Base.append!(sse::SSE, esses::Vector{ESSE}) = append!(sse.sequence, esses)
Base.append!(sse::SSE, esses::SSE) = append!(sse.sequence, esses.sequence)
Base.delete!(sse::SSE, i::Int) = deleteat!(sse.sequence, i)
Base.deleteat!(sse::SSE, i::Int) = deleteat!(sse.sequence, i)
Base.deleteat!(sse::SSE, inds) = deleteat!(sse.sequence, inds)
Base.splice!(sse::SSE, i::Int) = splice!(sse.sequence, i)
Base.splice!(sse::SSE, i::Int, esses::Vector{ESSE}) = splice!(sse.sequence, i, esses)
Base.keepat!(sse::SSE, inds) = keepat!(sse.sequence, inds)
Base.keepat!(sse::SSE, m::AbstractVector{Bool}) = keepat!(sse.sequence, m)

function Base.reverse!(sse::SSE)
    reverse!.(sse.sequence)
    return SSE(reverse!(sse.sequence))
end
function Base.reverse(sse::SSE)
    return reverse!(deepcopy(sse))
end

function Base.show(io::IO, sse::SSE)
    print(io, "StrongShiftEquivalence(")
    for esse in sse.sequence[1:end-1]
        print(io, esse[1] => esse[2])
        print(io, ", ")
    end
    print(io, sse.sequence[end][1] => sse.sequence[end][2])
    print(io, ")")
end

function Base.show(io::IO, ::MIME"text/plain", sse::SSE)
    println(io, "StrongShiftEquivalence with lag $(length(sse)):")
    d_h, d_w = displaysize(io)
    hs = Vector{Int}(undef, length(sse))
    ws = Vector{Int}(undef, length(sse))
    Rs_str = Vector{Vector{String}}(undef, length(sse))
    Ss_str = Vector{Vector{String}}(undef, length(sse))
    RSs_str = Vector{Vector{String}}(undef, length(sse))
    SRs_str = Vector{Vector{String}}(undef, length(sse))
    margin_w = 4
    margin_h = 4
    for (i, esse) in enumerate(sse)
        R, S = esse
        RS, SR = R * S, S * R
        R_str = split(matrix2str(R), "\n")
        S_str = split(matrix2str(S), "\n")
        RS_str = split(matrix2str(RS), "\n")
        SR_str = split(matrix2str(SR), "\n")
        hs[i] = max(length(R_str)+1, length(S_str)+1)
        if i == 1
            ws[i] = length(R_str[1]) + length(S_str[1]) + length(RS_str[1]) + length(SR_str[1]) + margin_w
        else
            ws[i] = length(R_str[1]) + length(S_str[1]) + length(SR_str[1]) + margin_w
        end
        Rs_str[i] = R_str
        Ss_str[i] = S_str
        RSs_str[i] = RS_str
        SRs_str[i] = SR_str
    end

    show_with_braille = false
    if maximum(hs) + margin_h > d_h || maximum(ws) + margin_w > d_w
        show_with_braille = true
    end

    if show_with_braille
        for (i, esse) in enumerate(sse)
            R, S = esse
            RS, SR = R * S, S * R
            R_str = split(matrix2braille(R), "\n")
            S_str = split(matrix2braille(S), "\n")
            RS_str = split(matrix2braille(RS), "\n")
            SR_str = split(matrix2braille(SR), "\n")

            hs[i] = max(length(R_str)+1, length(S_str)+1)
            if i == 1
                ws[i] = length(R_str[1]) + length(S_str[1]) + length(RS_str[1]) + length(SR_str[1]) + margin_w
            else
                ws[i] = length(R_str[1]) + length(S_str[1]) + length(SR_str[1]) + margin_w
            end

            Rs_str[i] = R_str
            Ss_str[i] = S_str
            RSs_str[i] = RS_str
            SRs_str[i] = SR_str
        end
        
        s = 1.0
        if maximum(hs) + margin_h > d_h || maximum(ws) + margin_w > d_w # + 4 is just for safety 
            s = min(d_h / (maximum(hs) + margin_h), d_w / (maximum(ws) + margin_w))
        end
        for (i, esse) in enumerate(sse)
            R, S = esse
            RS, SR = R * S, S * R
            RS_sw = floor(Int, length(RSs_str[i][1]) * s)
            RS_sh = floor(Int, length(RSs_str[i]) * s)
            SR_sw = floor(Int, length(SRs_str[i][1]) * s)
            SR_sh = floor(Int, length(SRs_str[i]) * s)

            Rs_str[i] = split(matrix2braille(R, RS_sh, nothing), "\n")
            Ss_str[i] = split(matrix2braille(S, SR_sh, nothing), "\n")
            RSs_str[i] = split(matrix2braille(RS, RS_sh, RS_sw), "\n")
            SRs_str[i] = split(matrix2braille(SR, SR_sh, SR_sw), "\n")
            hs[i] = max(length(Rs_str[i])+1, length(Ss_str[i])+1) #, length(RSs_str[i]), length(SRs_str[i]))
            if i > 1
                ws[i] = length(Rs_str[i][1]) + length(Ss_str[i][1]) + length(SRs_str[i][1]) + margin_w
            else
                ws[i] = length(Rs_str[i][1]) + length(Ss_str[i][1]) + length(RSs_str[i][1]) + length(SRs_str[i][1]) + margin_w
            end
        end
    end

    # construct string vectors for each esse
    local esse_strs = Vector{Vector{String}}(undef, length(sse))
    for (l, esse) in enumerate(sse)
        R_str = Rs_str[l]
        S_str = Ss_str[l]
        RS_str = RSs_str[l]
        SR_str = SRs_str[l]
        esse_str = Vector{String}(undef, hs[l])
        for i in 1:hs[l]-1
            esse_str[i] = ""
            if l == 1
                if i - (hs[l] - length(RS_str)) > 0
                    j = i - (hs[l] - length(RS_str))
                    esse_str[i] *= RS_str[j]
                else
                    esse_str[i] *=" "^(length(RS_str[1]))
                end
            end
            esse_str[i] *= " "
            if i - (hs[l] - length(R_str)) + 1 > 0
                j = i - (hs[l] - length(R_str)) + 1
                esse_str[i] *=R_str[j]
            else
                esse_str[i] *= " "^(length(R_str[1]))
            end
            esse_str[i] *= " "
            if i - (hs[l] - length(S_str)) + 1 > 0
                j = i - (hs[l] - length(S_str)) + 1
                esse_str[i] *= S_str[j]
            else
                esse_str[i] *= " "^(length(S_str[1]))
            end
            esse_str[i] *= "  "
            if i - (hs[l] - length(SR_str)) > 0
                j = i - (hs[l] - length(SR_str))
                esse_str[i] *= SR_str[j]
            else
                esse_str[i] *= " "^(length(SR_str[1]))
            end
        end
        esse_str[end] = l == 1 ? RS_str[end] : ""
        esse_str[end] *= "─"^(length(R_str[1]) + length(S_str[1]) + 3)
        esse_str[end] *= "▶︎" * SR_str[end]
        esse_strs[l] = esse_str
    end

    # print
    l = 0 # counter for drawn esses
    r = 1 # counter for drawn rows
    while l < length(sse)
        local draw_esses::UnitRange{Int}
        w = 0
        l += 1
        for i in l:length(sse)
            w += ws[i]
            if w + 4 > d_w # + 4 is just for safety
                draw_esses = i > l ? (l:i-1) : l:l
                l = last(draw_esses)
                break
            end
            if i == length(sse)
                draw_esses = l:length(sse)
                l = length(sse)
            end
        end

        max_h = maximum(hs[draw_esses])
        for i in reverse(2:max_h)
            for j in draw_esses
                j ∈ first(draw_esses) && r > 1 && print(io, "  ")
                if i <= length(esse_strs[j])
                    i_tmp = length(esse_strs[j]) - i + 1
                    print(io, esse_strs[j][i_tmp])
                else
                    print(io, " "^(length(esse_strs[j][1])))
                end
            end
            println(io)
        end
        (r > 1) && print(io, "⋯─")
        for j in draw_esses
            print(io, esse_strs[j][end])
        end
        l < length(sse) && println(io, "─⋯\n")

        r += max_h
    end
end

"""
    reduce!(sse::SSE)

Remove identity ESSE components wand return it.

# Examples
```jldoctest
f = extend(shift([2;;]), 1, 0)
2-BlockMap: SFT([2;;]) → SFT([2;;]) with memory=0:
(1, 1) => 1
(1, 2) => 2
(2, 1) => 1
(2, 2) => 2

sse = _kichens_StrongShiftEquivalence(f)[1]
StrongShiftEquivalence with lag 6:
                               ⎡1 0⎤            ⎡1 1 0 0⎤ ⎡1 0 0 0⎤
                               ⎢0 1⎥  ⎡1 1 0 0⎤ ⎢0 0 1 1⎥ ⎢0 1 0 0⎥  ⎡1 1 0 0⎤
          ⎡1⎤        ⎡1 1 0 0⎤ ⎢1 0⎥  ⎢0 0 1 1⎥ ⎢1 1 0 0⎥ ⎢0 0 1 0⎥  ⎢0 0 1 1⎥
    [1 1] ⎣1⎦  ⎡1 1⎤ ⎣0 0 1 1⎦ ⎣0 1⎦  ⎢1 1 0 0⎥ ⎣0 0 1 1⎦ ⎣0 0 0 1⎦  ⎢1 1 0 0⎥
[2]───────────▶︎⎣1 1⎦─────────────────▶︎⎣0 0 1 1⎦─────────────────────▶︎⎣0 0 1 1⎦─⋯

  ⎡1 1 0 0⎤ ⎡1 0 0 0⎤            ⎡1 0⎤
  ⎢0 0 1 1⎥ ⎢0 1 0 0⎥  ⎡1 1 0 0⎤ ⎢0 1⎥
  ⎢1 1 0 0⎥ ⎢0 0 1 0⎥  ⎢0 0 1 1⎥ ⎢1 0⎥ ⎡1 1 0 0⎤        ⎡1⎤
  ⎣0 0 1 1⎦ ⎣0 0 0 1⎦  ⎢1 1 0 0⎥ ⎣0 1⎦ ⎣0 0 1 1⎦  ⎡1 1⎤ ⎣1⎦ [1 1]
⋯─────────────────────▶︎⎣0 0 1 1⎦─────────────────▶︎⎣1 1⎦───────────▶︎[2]

reduce!(sse)
2-element Vector{ElementaryStrongShiftEquivalence}:
 ElementaryStrongShiftEquivalence([1 1 0 0; 0 0 1 1; 1 1 0 0; 0 0 1 1] => [1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 1])
 ElementaryStrongShiftEquivalence([1 1 0 0; 0 0 1 1; 1 1 0 0; 0 0 1 1] => [1 0 0 0; 0 1 0 0; 0 0 1 0; 0 0 0 1])

sse
StrongShiftEquivalence with lag 4:
                               ⎡1 0⎤            ⎡1 0⎤                                
                               ⎢0 1⎥  ⎡1 1 0 0⎤ ⎢0 1⎥                                
          ⎡1⎤        ⎡1 1 0 0⎤ ⎢1 0⎥  ⎢0 0 1 1⎥ ⎢1 0⎥ ⎡1 1 0 0⎤        ⎡1⎤           
    [1 1] ⎣1⎦  ⎡1 1⎤ ⎣0 0 1 1⎦ ⎣0 1⎦  ⎢1 1 0 0⎥ ⎣0 1⎦ ⎣0 0 1 1⎦  ⎡1 1⎤ ⎣1⎦ [1 1]     
[2]───────────▶︎⎣1 1⎦─────────────────▶︎⎣0 0 1 1⎦─────────────────▶︎⎣1 1⎦───────────▶︎[2]
```
"""
function reduce!(sse::SSE)
    removed = Vector{ESSE}(undef, 0)
    remove_inds = Vector{Int}(undef, 0)

    function is_identity(A::Matrix{Int})
        return A == oneunit(A)
    end

    if length(sse) <= 1
        return removed
    end

    for (i, (R, S)) in enumerate(sse)
        if size(R) != size(S)
            continue
        end
        if is_identity(R) || is_identity(S)
            push!(removed, sse[i])
            push!(remove_inds, i)
        end
    end

    deleteat!(sse, remove_inds)
    return removed
end

reduce(sse::SSE) = begin
    sse = deepcopy(sse)
    reduce!(sse)
    return sse
end

"""
    reducefirst!(sse::SSE)

Remove the first identity ESSE component and return it.
"""
function reducefirst!(sse::SSE)
    function is_identity(A::Matrix{Int})
        return A == oneunit(A)
    end
    if length(sse) <= 1
        return nothing
    end
    for (i, (R, S)) in enumerate(sse)
        if size(R) != size(S)
            continue
        end
        if is_identity(R) || is_identity(S)
            return popat!(sse, i)
        end
    end
    return nothing
end

function reducefirst(sse::SSE)
    sse = deepcopy(sse)
    reducefirst!(sse)
    return sse
end

"""
    reducelast!(sse::SSE)

Remove the last identity ESSE component and return it.
"""
function reducelast!(sse::SSE)
    function is_identity(A::Matrix{Int})
        return A == oneunit(A)
    end
    if length(sse) <= 1
        return nothing
    end
    for (i, (R, S)) in enumerate(reverse(sse))
        if size(R) != size(S)
            continue
        end
        if is_identity(R) || is_identity(S)
            return popat!(sse, length(sse) - i + 1)
        end
    end
    return nothing
end

function reducelast(sse::SSE)
    sse = deepcopy(sse)
    reducelast!(sse)
    return sse
end

function _vertex_identity(A::Matrix{Int})
    blocks = Dict{NTuple{1, Int}, Int}()
    n = size(A)[1]

    for i in 1:n
        blocks[(i,)] = i
    end

    block_map = BlockMap(A, A, blocks, 0)
    return block_map
end

function _edge_identity(A::Matrix{Int})
    blocks = Dict{NTuple{1, Int}, Int}()

    for i in 1:sum(A)
        blocks[(i,)] = i
    end

    block_map = BlockMap(A, A, blocks, 0)
    return block_map
end

"""
    identity(A::Matrix{Int})

Return the identity block map of the given matrix.

# Examples
```jldoctest
```
"""
function identity(A::Matrix{Int})
    if size(A)[1] != size(A)[2]
        throw(ArgumentError("matrix must be a square"))
    end

    if all(0 .<= A .<= 1)
        return _vertex_identity(A)
    else
        return _edge_identity(A)
    end
end

function Base.oneunit(block_map::BlockMap)
    if size(block_map.start)[1] != size(block_map.target)[1]
        throw(ArgumentError("The start and target matrix must have the same size"))
    end
    return identity(block_map.start)
end

"""
    shift(A::Matrix{Int})

Return the shift block map of the given matrix. 
If A contains only 0 and 1, the function returns the block map on the vertex shift. 
Otherwise, the function returns the block map on the edge shift.

# Examples
```jldoctest
```
"""
function shift(A::Matrix{Int})
    if size(A)[1] != size(A)[2]
        throw(ArgumentError("matrix must be a square"))
    end

    if all(0 .<= A .<= 1)
        return _vertex_shift(A)
    else
        return _edge_shift(A)
    end
end

function _vertex_shift(A::Matrix{Int})
    blocks = Dict{Tuple{Int}, Int}()
    n = size(A)[1]

    for i in 1:n
        blocks[(i,)] = i
    end

    return BlockMap(A, A, blocks, -1)
end

function _edge_shift(A::Matrix{Int})
    blocks = Dict{Tuple{Int}, Int}()
    s = sum(A)

    for i in 1:s
        blocks[(i,)] = i
    end

    return BlockMap(A, A, blocks, -1)
end

"""
    marker_automorphism(markers::Tuple{Vararg{Array{Union{Int, Char}, 1}}})

Return the compound marker automorphism on the full two shift.
"""
function marker_endomorphism(markers::Tuple{Vararg{Array{Union{Int, Char}, 1}}})
    # markers は数字または'*'からなる配列のタプルである
    # 例えば ([1, '*', 2, 1], [2, '*', 1, 1]) など
    # 各マーカーには必ずただ1つの '*' が含まれている

    # throw error if the asterisk is not included just once in each marker
    for marker in markers
        if length(findall(x -> x == '*', marker)) != 1
            throw(ArgumentError("Each marker must have just one *"))
        end
    end

    # record position of *
    pos = Vector{Int}()
    for marker in markers
        push!(pos, findfirst(x -> x == '*', marker))
    end
    max_pos = maximum(pos)
    max_before_len = maximum(pos) - 1

    # record the length of each marker
    len = Vector{Int}()
    for marker in markers
        push!(len, length(marker))
    end
    max_after_len = maximum(len .- pos)

    # padding with '.' based on the position of *
    for (i, marker) in enumerate(markers)
        pushfirst!(marker, fill('.', max_before_len - pos[i] + 1)...)
        push!(marker, fill('.', max_after_len - (len[i] - pos[i]))...)
    end

    # constuct the blocks for the compound marker automorphism
    N = maximum(length.(markers))
    blocks = Dict{NTuple{N, Int}, Int}()
    for w in EdgeWalks([2;;], maximum(length.(markers)))
        y = w[max_pos]
        for marker in markers
            is_hit = true
            for (i, c) in enumerate(marker)
                if !(c ∈ ".*") && c != w[i]
                    is_hit = false
                    break
                end
            end
            if is_hit
                # flipping 1 and 2 at the position of *
                y = -w[max_pos] + 3
                break
            end
        end
        push!(blocks, w => y)
    end

    return shorten(BlockMap([2;;], [2;;], blocks, max_pos-1))
end


"""
    marker_automorphism(markers::Tuple{Vararg{String}})

Return the compound marker automorphism on the full two shift.

# Examples
```jldoctest
julia> marker_automorphism("1 1 * 2", "2 * 2 1")
BlockMap([2;;], [2;;], Dict{Tuple{Vararg{Int64}}, Int64}((1, 1, 2, 2, 2) => 1, (1, 1, 2, 1, 2) => 2, (1, 2, 2, 2, 1) => 1, (1, 2, 2, 1, 1) => 2, (2, 1, 1, 2, 1) => 1, (2, 1, 1, 1, 1) => 1, (1, 2, 2, 2, 2) => 2, (1, 2, 2, 1, 2) => 2, (2, 1, 1, 2, 2) => 1, (2, 1, 1, 1, 2) => 1…), 2)
"""
function marker_endomorphism(markers::Tuple{Vararg{String}})
    splited = Array{Array{Union{Int, Char}, 1}, 1}()
    for marker in markers
        tmp = split(marker, " ")
        tmp = [x == "*" ? '*' : parse(Int, x) for x in tmp]
        push!(splited, tmp)
    end
    return marker_endomorphism(tuple(splited...))
end

marker_endomorphism(markers::Vector{String}) = marker_endomorphism(tuple(markers...))
marker_endomorphism(markers::String...) = marker_endomorphism(markers)

"transform edge shift matrix to vertex one and compute ESSE of the transformation"
function _vertex(A::Matrix{Int})
    if size(A)[1] != size(A)[2]
        throw(ArgumentError("matrix must be a square"))
    end

    ind = cumsum([0, A...])
    n = size(A)[1]

    R = zeros(Int, n, sum(A))
    S = zeros(Int, sum(A), n)

    for i in 0:n-1
        S[ind[i*n+1]+1:ind[(i+1)*n+1], i+1] .= 1
        for ii in 0:n-1
            R[i+1, ind[(i+1) + ii * n]+1:ind[(i+1) + ii * n + 1]] .= 1
        end
    end
    
    return S * R, ESSE(R => S)
end

"calculate the higher order representation of the given block map"
function _higher_order_representation(block_map::BlockMap, x_order::Int, y_order::Int)
    on_s, _ = _is_on_vertex_shift(block_map)
    if on_s
        return _hi_ord_rep_vertex(block_map, x_order, y_order)
    else
        return _hi_ord_rep_edge(block_map, x_order, y_order)
    end
end

function _hi_ord_rep_vertex(block_map::BlockMap, x_order::Int, y_order::Int)
    w_l = max(x_order, length(block_map) + y_order - 1)
    L = Set{Pair{NTuple{x_order, Int}, NTuple{y_order, Int}}}()

    for w in VertexWalks(block_map.start, w_l)
        push!(L, w[1:x_order] => block_map(w)[1:y_order])
    end

    return collect(L) |> unique
end

function _hi_ord_rep_edge(block_map::BlockMap, x_order::Int, y_order::Int)
    w_l = max(x_order, length(block_map) + y_order - 1)
    L = Set{Pair{NTuple{x_order, Int}, NTuple{y_order, Int}}}()

    for w in EdgeWalks(block_map.start, w_l)
        push!(L, w[1:x_order] => block_map(w)[1:y_order])
    end

    return collect(L) |> unique
end

function (g::BlockMap)(f::BlockMap)
    if f.target != g.start
        throw(ArgumentError("block maps must be composable"))
    end

    l1 = length(f)
    l2 = length(g)
    blocks = Dict{NTuple{l1 + l2 - 1, Int}, Int}()
    if _is_on_vertex_shift(f)[1]
        for w in VertexWalks(f.start, l1 + l2 - 1)
            blocks[w] = g(f(w))[1]
        end
    else
        for w in EdgeWalks(f.start, l1 + l2 - 1)
            blocks[w] = g(f(w))[1]
        end
    end

    return shorten(BlockMap(f.start, g.target, blocks, f.memory + g.memory))
end

"""
    inv(block_map::BlockMap)

Return the inverse block map of the given block map.

# Examples
```jldoctest
```
"""
function Base.inv(block_map::BlockMap)
    block_map = shorten(block_map)
    l = length(block_map)
    
    for i in 1:l+5
        L = _higher_order_representation(block_map, max(l, i), i)
        for j in 1:i
            blocks = Dict{typeof(L[1][2]), Int}()
            find_inv = true
            for (x, y) in L
                if !haskey(blocks, y)
                    blocks[y] = x[j]
                elseif blocks[y] != x[j]
                    # not injective case
                    find_inv = false
                    break
                end
            end
            if find_inv
                try
                    return shorten(BlockMap(block_map.start, block_map.target, blocks, -block_map.memory + j - 1))
                catch e
                    # not surjective case
                    continue
                end
            end
        end
    end

    throw(ArgumentError("we could not find the inverse"))
end

"""
    (f::BlockMap) ^ (n::Int)

Compute the n-th power of the given block map.
`n` is able to be negative if the map is invertible.
"""
function (Base.:^)(block_map::BlockMap, n::Int)
    if block_map.start != block_map.target
        throw(ArgumentError("map is not self composable"))
    end

    if n >= 1
        for i in 1:n-1
            block_map = block_map(block_map)
        end
        return block_map
    elseif n == 0
        return identity(block_map.start)
    else # n < 0
        block_map = inv(block_map)
        for i in 1:(-n-1)
            block_map = block_map(block_map)
        end 
        return inv(block_map)
    end
end

"""
    extend(bloc_map::BlockMap, pre::Int, post::Int) -> ::BlockMap

Extend the defining blocks of given block map.

This function return a block map which has definition blocks extended by `pre` in forward and `post` in backward.

If `pre` or `post` is negative, this function runs `shorten` instead.

# Examples
"""
function extend(block_map::BlockMap, pre::Int = 0, post::Int = 0)
    nonneg_pre = max(0, pre)
    nonneg_post = max(0, post)

    blockWalks = _is_on_vertex_shift(block_map)[1] ? VertexWalks : EdgeWalks

    blocks = Dict{NTuple{length(block_map) + nonneg_pre + nonneg_post, Int}, Int}()
    for w in blockWalks(block_map.start, nonneg_pre + length(block_map) + nonneg_post)
        blocks[w] = block_map(w[nonneg_pre+1:end-nonneg_post])[1]
    end

    if pre < 0 || post < 0
        neg_pre = max(0, -pre)
        neg_post = max(0, -post)
        return shorten(BlockMap(block_map.start, block_map.target, blocks, block_map.memory + pre), neg_pre, neg_post)
    end
    return BlockMap(block_map.start, block_map.target, blocks, block_map.memory + pre)
end

"""
    shortenable(block_map::BlockMap) -> ::Tuple{Int, Int}

Return the number of blocks which can be shortened in forward and backward.

# Examples
```jldoctest
f = extend(shift([2;;]), 1, 2)
4-BlockMap: SFT([2;;]) → SFT([2;;]) with memory=0:
(1, 1, 1, 1) => 1
(1, 1, 1, 2) => 1
(1, 1, 2, 1) => 1
(1, 1, 2, 2) => 1
        ⋮
(2, 2, 1, 1) => 2
(2, 2, 1, 2) => 2
(2, 2, 2, 1) => 2
(2, 2, 2, 2) => 2

julia> shortenable(f)
(1, 2)
"""
function shortenable(block_map::BlockMap)
    blocks = block_map.blocks
    blockWalks = _is_on_vertex_shift(block_map)[1] ? VertexWalks : EdgeWalks

    pre = length(block_map) - 1
    for i in 1:length(block_map) - 1, w in blockWalks(block_map.start, length(block_map) - i)
        idx = findall(x -> x[1][i+1:end] == w, collect(blocks))
        ys = [y for (_, y) in collect(blocks)[idx]]
        if length(Set(ys)) > 1
            pre = i - 1
            break
        end
    end

    post = length(block_map) - pre - 1
    for i in 1:length(block_map) - pre - 1, w in blockWalks(block_map.start, length(block_map) - i)
        idx = findall(x -> x[1][1:end-i] == w, collect(blocks))
        ys = [y for (_, y) in collect(blocks)[idx]]
        if length(Set(ys)) > 1
            post = i - 1
            break
        end
    end

    return pre, post
end

"""
    shorten(block_map::BlockMap)

Return the block map which is same as the given block map but the length of defining blocks are the shortest.

# Examples
```jldoctest
f = extend(shift([2;;]), 1, 2)
4-BlockMap: SFT([2;;]) → SFT([2;;]) with memory=0:
(1, 1, 1, 1) => 1
(1, 1, 1, 2) => 1
(1, 1, 2, 1) => 1
(1, 1, 2, 2) => 1
        ⋮
(2, 2, 1, 1) => 2
(2, 2, 1, 2) => 2
(2, 2, 2, 1) => 2
(2, 2, 2, 2) => 2

shorten(f)
1-BlockMap: SFT([2;;]) → SFT([2;;]) with memory=-1:
(1,) => 1
(2,) => 2
```
"""
function shorten(block_map::BlockMap)
    pre, post = shortenable(block_map)

    blocks = Dict{NTuple{length(block_map) - pre - post, Int}, Int}()
    for (x, y) in block_map.blocks
        blocks[x[pre+1:end-post]] = y
    end

    return BlockMap(block_map.start, block_map.target, blocks, block_map.memory - pre)
end

"""
    shorten(block_map::BlockMap, pre::Int, post::Int)

Return the block map which is same as the given block map but the length of defining blocks are shorten specific length in front and back.

`pre` and `post` must be less than `shortenable(block_map)`.

If `pre` or `post` is negative, this function runs `extend` instead.

# Examples
```jldoctest
f = extend(shift([2;;]), 1, 2)
4-BlockMap: SFT([2;;]) → SFT([2;;]) with memory=0:
(1, 1, 1, 1) => 1
(1, 1, 1, 2) => 1
(1, 1, 2, 1) => 1
(1, 1, 2, 2) => 1
        ⋮
(2, 2, 1, 1) => 2
(2, 2, 1, 2) => 2
(2, 2, 2, 1) => 2
(2, 2, 2, 2) => 2

shorten(f, 1, 1)
2-BlockMap: SFT([2;;]) → SFT([2;;]) with memory=-1:
(1, 1) => 1
(1, 2) => 1
(2, 1) => 2
(2, 2) => 2
````
"""
function shorten(block_map::BlockMap, pre::Int, post::Int)
    nonneg_pre = max(0, pre)
    nonneg_post = max(0, post)

    max_pp = shortenable(block_map)

    if max_pp[1] < pre
        @warn "pre is exceed the maximum shortenable length"
    end
    if max_pp[2] < post
        @warn "post is exceed the maximum shortenable length"
    end

    blocks = Dict{NTuple{length(block_map) - nonneg_pre - nonneg_post, Int}, Int}()
    for (x, y) in block_map.blocks
        blocks[x[nonneg_pre+1:end-nonneg_post]] = y
    end

    if pre < 0 || post < 0
        neg_pre = max(0, -pre)
        neg_post = max(0, -post)
        return extend(BlockMap(block_map.start, block_map.target, blocks, block_map.memory - pre), neg_pre, neg_post)
    end
    return BlockMap(block_map.start, block_map.target, blocks, block_map.memory - pre)
end

"""
    Matrix(block_map::BlockMap)

Return the matrix representation of the given block map.

# Examples
```jldoctest
```
"""
function Base.Matrix(block_map::BlockMap)
    l = length(block_map)
    L = _higher_order_representation(block_map, l, l)
    words = [x for (x, _) in L]
    unique!(words)

    M = zeros(Int, length(words), length(words))
    for (x, y) in L
        M[findfirst(u -> u == x, words), findfirst(v -> v == y, words)] = 1
    end

    return M
end

"""
    StrongShiftEquivalence(::BlockMap) -> ::Tuple{SSE{Int}, Int}

Compute the strong shift equivalence over ℤ₊ (SSE-ℤ₊) of the given block map and how many times the δ(A)⁻¹ must be applied to get the dimension representation of the block map.

You can also use `SSE(block_map::BlockMap)` instead of `StrongShiftEquivalence(block_map::BlockMap)`.

# Examples
```jldoctest
```
"""
function StrongShiftEquivalence(block_map::BlockMap; alg::Symbol = :Kichens)
    if alg == :Kichens
        return _kichens_StrongShiftEquivalence(block_map)
    # #TODO: implement below
    # elseif alg == :LindMarcus
    #     return _lindmarcus_StrongShiftEquivalence(block_map)
    else
        throw(ArgumentError("algorithm must be :Kichens"))
    end
end

function _kichens_StrongShiftEquivalence(block_map::BlockMap)
    # TODO: fix this function
    # TODO: modify for the case that memory is lager than length(block_map)
    # TODO: implement both of Kitchen's and Lind Marcus' version
    # TODO: in this scope, make function that construct _SR
    
    # Step 0 - prepare the information
    inv_map = inv(block_map)

    A = block_map.start
    B = block_map.target

    n = size(A)[1]
    m = size(B)[1]
    l = length(block_map)
    l_inv = length(inv_map)
    memory = block_map.memory
    memory_inv = inv_map.memory

    list_Sf = Matrix{Int}[]
    list_Rf = Matrix{Int}[]
    list_Sp = Matrix{Int}[]
    list_Rp = Matrix{Int}[]

    # `num_R_not_subdiv` count the number of times which matrix R is not a subdivision matrix
    num_R_not_subdiv = 0

    Revf = Matrix(I, n, n)
    Sevf = Matrix(I, n, n)
    Revp = Matrix(I, m, m)
    Sevp = Matrix(I, m, m)

    # convert edge shift to vertex shift
    is_vertex = _is_on_vertex_shift(block_map)
    if !is_vertex[1]
        A, (Revf, Sevf) = _vertex(A)
        n = size(A)[1]
    end
    if !is_vertex[2]
        B, (Revp, Sevp) = _vertex(B)
        m = size(B)[1]
        num_R_not_subdiv += 1
    end

    # compute shift number 
    a = memory + 1 - l
    b = memory
    b_inv = -memory_inv
    a_inv = -(memory_inv + 1 - l_inv)
    # check if intervals [a, b] and [a_inv, b_inv] are intersected
    if intersect(Set(a:b), Set(a_inv:b_inv)) != Set()
        num_applied_shift = minimum(intersect(Set(a:b), Set(a_inv:b_inv)))
    else
        k = 0
        if b_inv < a
            k = round(Int, (b_inv + a) / 2)
        else
            k = round(Int, (b + a_inv) / 2)
        end
        forward_length(k) = max(l, memory - k, -memory + k + l)
        backward_length(k) = max(l_inv, memory_inv + k, -memory_inv - k + l_inv)
        num_applied_shift = argmin(
            x -> length(VertexWalks(A, forward_length(x))) + length(VertexWalks(B, backward_length(x))), 
            [k, a, b, a_inv, b_inv]
        )
    end
    x_zero_index = memory + 1 - num_applied_shift
    y_zero_index = memory_inv + 1 + num_applied_shift

    # extend the block map if x_zero_index (y_zero_index) is out of the range [1, l] ([1, l_inv])
    _pre = x_zero_index < 1 ? -x_zero_index + 1 : 0
    _post = x_zero_index > l ? x_zero_index - l : 0
    block_map = extend(block_map, _pre, _post)

    _inv_pre = y_zero_index < 1 ? -y_zero_index + 1 : 0
    _inv_post = y_zero_index > l_inv ? y_zero_index - l_inv : 0
    inv_map = extend(inv_map, _inv_pre, _inv_post)

    # sort by order that makes subdivision matrix like a diagonal matrix
    full_blocks = Set{Pair{NTuple{l, Int}, NTuple{l_inv, Int}}}()
    for w in VertexWalks(A, l + l_inv - 1)
        push!(full_blocks, w[y_zero_index:y_zero_index + l - 1] => block_map(w))
    end
    full_blocks = collect(full_blocks)
    sort!(full_blocks, by = xy -> [
        reverse(xy[1][1:x_zero_index]) ;
        xy[1][x_zero_index+1:end] ;
        reverse(xy[2][1:y_zero_index]) ;
        xy[2][y_zero_index+1:end]
    ])
    blocks = [x => y[y_zero_index] for (x, y) in full_blocks] |> unique

    _SR = B

    # step 1 - split A to higher block representation in x
    # step 1.1 - extend to past direction
    L = blocks
    if x_zero_index > 1
        # construct x_zero_index-length symbols
        new_L = [x[1:x_zero_index] for (x, y) in blocks] |> unique
        # split 方向
        # | S : subdivision matrix(縦長行列)
        # | R : 横長行列

        # construct SR
        _SR = zeros(Int, length(new_L), length(new_L))
        for (i, x1) in enumerate(new_L), (j, x2) in enumerate(new_L)
            if x1[2:end] == x2[1:end-1]# && y1[2:end] == y2[1:end-1]
                _SR[i, j] = 1
            end
        end

        # amalgamate the symbols step by step
        tmp_list_Sf = Matrix{Int}[]
        tmp_list_Rf = Matrix{Int}[]
        for _ in 2:x_zero_index
            L = [x[2:end] for x in new_L] |> unique
            _split = Vector{Vector{Int}}(undef, length(L))
            fill!(_split, [])
            for (i, x) in enumerate(new_L)
                j = findfirst(_x -> _x == x[2:end], L)
                tmp = copy(_split[j])
                push!(tmp, i)
                _split[j] = tmp
            end

            S = zeros(Int, length(new_L), length(L))
            R = zeros(Int, length(L), length(new_L))

            # construct S
            for i in axes(S, 2)
                S[_split[i], i] .= 1
            end

            # construct R
            for i in axes(R, 1)
                R[i, :] = _SR[_split[i][1], :]
            end

            _SR = R * S

            pushfirst!(tmp_list_Sf, S)
            pushfirst!(tmp_list_Rf, R)

            # update symbol set
            new_L = L
        end
        append!(list_Sf, tmp_list_Sf)
        append!(list_Rf, tmp_list_Rf)
        num_R_not_subdiv += length(tmp_list_Rf)
    else
        L = [(x, ) for x in 1:n]
    end

    # step1.2 extend to future direction
    # S : 縦長行列
    # R : subdivision matrix(横長行列)
    L = [x[1:x_zero_index] for (x, y) in blocks] |> unique
    if x_zero_index < l
        new_L = [x for (x, y) in blocks]

        # construct SR
        _SR = zeros(Int, length(new_L), length(new_L))
        for (i, x1) in enumerate(new_L), (j, x2) in enumerate(new_L)
            if x1[2:end] == x2[1:end-1]# && y1[2:end] == y2[1:end-1]
                _SR[i, j] = 1
            end
        end

        # amalgamate the symbols step by step
        tmp_list_Rf = Matrix{Int}[]
        tmp_list_Sf = Matrix{Int}[]
        for _ in x_zero_index+1:l
            L = [x[1:end-1] for x in new_L] |> unique
            _split = Vector{Vector{Int}}(undef, length(L)) # pair of index of L and split symbols
            fill!(_split, [])
            for (i, x) in enumerate(new_L)
                j = findfirst(_x -> _x == x[1:end-1], L)
                tmp = copy(_split[j])
                push!(tmp, i)
                _split[j] = tmp
            end

            S = Matrix{Int}(undef, length(new_L), length(L))
            R = zeros(Int, length(L), length(new_L))

            # construct R
            for i in axes(R, 1)
                R[i, _split[i]] .= 1
            end

            # construct S
            for i in axes(S, 2)
                S[:, i] = _SR[:, _split[i][1]]
            end

            _SR = R * S

            pushfirst!(tmp_list_Rf, R)
            pushfirst!(tmp_list_Sf, S)

            new_L = L
        end

        append!(list_Sf, tmp_list_Sf)
        append!(list_Rf, tmp_list_Rf)
    else
        # sort block_map.blocks that will be compatible with L
        # perhaps, this is not necessary because blocks is already sorted and L, new_L are constructed from blocks
        #=
        new_L = Vector{Pair{NTuple{l, Int}, Int}}(undef, length(blocks))
        for (i, x) in enumerate(L)
            j = findfirst(_xy -> _xy[1] == x, blocks)
            new_L[i] = blocks[j]
        end
        =#
    end

    L = blocks

    # step 2 - split image symbols y
    # step 2.1 - split in past direction
    # S : subdivision matrix (縦長行列)
    # R : 横長行列
    if y_zero_index > 1
        new_L = [x => y[1:y_zero_index] for (x, y) in full_blocks] |> unique

        # construct SR
        _SR = zeros(Int, length(new_L), length(new_L))
        for (i, (x1, y1)) in enumerate(new_L), (j, (x2, y2)) in enumerate(new_L)
            if x1[2:end] == x2[1:end-1] && y1[2:end] == y2[1:end-1]
                _SR[i, j] = 1
            end
        end

        # amalgamate the symbols step by step
        tmp_list_Rf = Matrix{Int}[]
        tmp_list_Sf = Matrix{Int}[]
        for _ in 2:y_zero_index
            L = [(x, y[2:end]) for (x, y) in new_L] |> unique
            _split = Vector{Vector{Int}}(undef, length(L)) # pair of index of L0 and split symbols
            fill!(_split, [])
            for (i, (x, y)) in enumerate(new_L)
                j = findfirst(xy -> xy[1] == x && xy[2] == y[2:end], L)
                tmp = copy(_split[j])
                push!(tmp, i)
                _split[j] = tmp
            end

            S = zeros(Int, length(new_L), length(L))
            R = Matrix{Int}(undef, length(L), length(new_L))
            # construct S
            for i in axes(S, 2)
                S[_split[i], i] .= 1
            end
        
            # construct R
            for i in axes(R, 1)
                R[i, :] = _SR[_split[i][1], :]
            end

            _SR = R * S

            pushfirst!(tmp_list_Sf, S)
            pushfirst!(tmp_list_Rf, R)

            # update symbol set
            new_L = L
        end
        append!(list_Sf, tmp_list_Sf)
        append!(list_Rf, tmp_list_Rf)
        num_R_not_subdiv += length(tmp_list_Rf)
    else
        nothing # to do
    end

    # step 2.2 - split in future direction
    # S : 縦長行列
    # R : subdivision matrix (横長行列)
    if y_zero_index < l_inv
        new_L = full_blocks
        # construct SR
        _SR = zeros(Int, length(new_L), length(new_L))
        for (i, (x1, y1)) in enumerate(new_L), (j, (x2, y2)) in enumerate(new_L)
            if x1[2:end] == x2[1:end-1] && y1[2:end] == y2[1:end-1]
                _SR[i, j] = 1
            end
        end

        # amalgamate the symbols step by step
        tmp_list_Rf = Matrix{Int}[]
        tmp_list_Sf = Matrix{Int}[]
        for _ in y_zero_index+1:l_inv
            L = [(x, y[1:end-1]) for (x, y) in new_L] |> unique
            _split = Vector{Vector{Int}}(undef, length(L)) # pair of index of L0 and split symbols
            fill!(_split, [])
            for (i, (x, y)) in enumerate(new_L)
                j = findfirst(xy -> xy[1] == x && xy[2] == y[1:end-1], L)
                tmp = copy(_split[j])
                push!(tmp, i)
                _split[j] = tmp
            end

            S = Matrix{Int}(undef, length(new_L), length(L))
            R = zeros(Int, length(L), length(new_L))

            # construct R
            for i in axes(R, 1)
                R[i, _split[i]] .= 1
            end
            
            # construct S
            for i in axes(S, 2)
                S[:, i] = _SR[:, _split[i][1]]
            end

            _SR = R * S

            pushfirst!(tmp_list_Sf, S)
            pushfirst!(tmp_list_Rf, R)

            # update symbol set
            new_L = L
        end
        append!(list_Sf, tmp_list_Sf)
        append!(list_Rf, tmp_list_Rf)
    else
        nothing # to do
    end

    # to be globally used
    !isempty(list_Sf) && (_SR = list_Sf[end] * list_Rf[end])

    # step 3 - amalgamate x symbols
    # step 3.1 - amalgamate in future direction
    # S : subdivision matrix (横長行列)
    # R : 縦長行列
    if x_zero_index < l
        new_L = full_blocks
        for _ in x_zero_index+1:l
            L = [x[1:end-1] => y for (x, y) in new_L] |> unique
            _split = Vector{Vector{Int}}(undef, length(L)) # pair of index of L and split symbols
            fill!(_split, [])
            for (i, (x, y)) in enumerate(new_L)
                j = findfirst(_xy -> _xy[1] == x[1:end-1] && _xy[2] == y, L)
                tmp = copy(_split[j])
                push!(tmp, i)
                _split[j] = tmp
            end

            S = zeros(Int, length(L), length(new_L))
            R = Matrix{Int}(undef, length(new_L), length(L))

            # construct S
            for i in axes(S, 1)
                S[i, _split[i]] .= 1
            end
            # construct R
            for i in axes(R, 2)
                R[:, i] = _SR[:, _split[i][1]]
            end

            _SR = S * R

            push!(list_Sf, S)
            push!(list_Rf, R)
            num_R_not_subdiv += 1

            new_L = L
        end
    else
        new_L = [x[1:x_zero_index] => y for (x, y) in full_blocks] |> unique
    end

    # step 3.2 - amalgamate in past direction
    # S : 横長行列
    # R : subdivision matrix (縦長行列)
    if x_zero_index > 1
        # new_L = [x[1:x_zero_index] => y for (x, y) in full_blocks] |> unique
        for _ in 2:x_zero_index
            L = [x[2:end] => y for (x, y) in new_L] |> unique
            _split = Vector{Vector{Int}}(undef, length(L))
            fill!(_split, [])
            for (i, (x, y)) in enumerate(new_L)
                j = findfirst(_xy -> _xy[1] == x[2:end] && _xy[2] == y, L)
                tmp = copy(_split[j])
                push!(tmp, i)
                _split[j] = tmp
            end

            S = Matrix{Int}(undef, length(L), length(new_L))
            R = zeros(Int, length(new_L), length(L))

            # construct R
            for i in axes(R, 2)
                R[_split[i], i] .= 1
            end
            # construct S
            for i in axes(S, 1)
                S[i, :] = _SR[_split[i][1], :]
            end

            _SR = S * R

            push!(list_Sf, S)
            push!(list_Rf, R)

            new_L = L
        end
    else
        nothing # to do
    end

    L = [y for (x, y) in L]

    # step 4 - amalgamate y symbols
    # step 4.1 - amalgamate in future direction (and sort by y[y_zero_index])
    # S : subdivision matrix (横長行列)
    # R : 縦長行列
    if y_zero_index < l_inv
        new_L = L
        for _ in y_zero_index+1:l_inv
            L = [y[1:end-1] for y in new_L] |> unique
            _split = Vector{Vector{Int}}(undef, length(L))
            fill!(_split, [])
            for (i, y) in enumerate(new_L)
                j = findfirst(_y -> _y == y[1:end-1], L)
                tmp = copy(_split[j])
                push!(tmp, i)
                _split[j] = tmp
            end

            S = zeros(Int, length(L), length(new_L))
            R = Matrix{Int}(undef, length(new_L), length(L))

            # construct S
            for i in axes(S, 1)
                S[i, _split[i]] .= 1
            end
            # construct R
            for i in axes(R, 2)
                R[:, i] = _SR[:, _split[i][1]]
            end

            _SR = S * R

            push!(list_Sf, S)
            push!(list_Rf, R)
            num_R_not_subdiv += 1

            new_L = L
        end
    else
        nothing # to do
    end

    # step 4.2 - amalgamate in past direction (and sort by y[y_zero_index])
    # S : 横長行列
    # R : subdivision matrix (縦長行列)
    if y_zero_index > 1
        new_L = L
        for _ in 2:y_zero_index
            L = [y[2:end] for y in new_L] |> unique
            _split = Vector{Vector{Int}}(undef, length(L))
            fill!(_split, [])
            for (i, y) in enumerate(new_L)
                j = findfirst(_y -> _y == y[2:end], L)
                tmp = copy(_split[j])
                push!(tmp, i)
                _split[j] = tmp
            end

            S = Matrix{Int}(undef, length(L), length(new_L))
            R = zeros(Int, length(new_L), length(L))

            # construct R
            for i in axes(R, 2)
                R[_split[i], i] .= 1
            end
            # construct S
            for i in axes(S, 1)
                S[i, :] = _SR[_split[i][1], :]
            end

            _SR = S * R

            push!(list_Sf, S)
            push!(list_Rf, R)

            new_L = L
        end
    else
        nothing # to do
    end

    # step 4.3 - compute permutation matrix
    P = zeros(Int, n, n)
    # skip if the permutation is identity
    _skip_flag = true
    for (i, y) in enumerate(L)
        P[i, y[1]] = 1
        i != y[1] && (_skip_flag = false)
    end
    !_skip_flag && push!(list_Rf, P)
    !_skip_flag && push!(list_Sf, B * P')

    # step 5 - calculate SSE that represent the block map
    list_esse = Array{ElementaryStrongShiftEquivalence{Int}, 1}(undef, 0)
    if !is_vertex[1]
        push!(list_esse, ESSE(Revf => Sevf))
    end
    for (R, S) in zip(list_Rf, list_Sf)
        push!(list_esse, ESSE(R => S))
    end
    for (S, R) in zip(list_Rp, list_Sp)
        push!(list_esse, ESSE(S => R))
    end
    if !is_vertex[2]
        push!(list_esse, ESSE(Sevp => Revp))
    end

    if length(list_esse) <= 0
        push!(list_esse, ESSE(Matrix{Int}(I, n, n) => B))
    end
    sse = StrongShiftEquivalence(list_esse)

    # (SSE, how many times the δ(A)⁻¹ must be applied)
    # `num_R_not_subdiv` is the number of times which matrix R is not a subdivision matrix
    return sse, num_R_not_subdiv + num_applied_shift
end

"""
    dimension_representation(SSE::StrongShiftEquivalence, [shift_inv::Int])

Compute the dimension representation of the given strong shift equivalence.
"""
function dimension_representation(sse::StrongShiftEquivalence, shift_inv::Int=0)
    A = sse[1][1] * sse[1][2]
    B = sse[end][2] * sse[end][1]
    M = oneunit(A)
    for esse in sse
        M *= esse[1]
    end
    if shift_inv > 0
        return M *= pinv(B)^shift_inv
    else
        return M * B^(-shift_inv)
    end
end

dimension_representation(x::Tuple{SSE, Int}) = dimension_representation(x[1], x[2])
dimension_representation(block_map::BlockMap) = dimension_representation(_kichens_StrongShiftEquivalence(block_map))

"""
    orbitsign_number(f::BlockMap, n)

Compute the orbit sign number of the given block map.

The orbit sign number is the sign of the block map as permutation of periodic orbits.
`os(f, n)` can be use for the value to be in ℤ/2ℤ.
"""
function orbitsign_number(f::BlockMap, n::Int)
    Os, Ot = shift_orbits(f, n, true)
    pdata = Vector{Int}(undef, length(Os))
    for (i, o) in enumerate(Os)
        # extend periodic sequence for the sequence f(o[1]) to be length n sequence
        in_seq = repeat([o[1]...], cld(length(f) - 1, n) + 1)[1:n+length(f)-1]
        pdata[i] = findfirst(x -> f(in_seq) in x, Ot)
    end
    return sign(Permutation(pdata))
end

"`0-1` version of `orbitsign_number`"
os(f, n) = 1 == orbitsign_number(f, n) ? 0 : 1

"""
    gyration_number(f::BlockMap, n)

Compute the gyration number of the given block map.
`gy(f, n)` also can be used.
"""
function gyration_number(f::BlockMap, n::Int)
    Os, Ot = shift_orbits(f, n, true)
    gyration = 0
    for o in Os
        in_seq = repeat([o[1]...], cld(length(f) - 1, n) + 1)[1:n+length(f)-1]
        out_arr = [f(in_seq)...]
        out_seq = tuple(out_arr...)
        o_ind = findfirst(x -> out_seq in x, Ot)
        for i in 0:n-1
            if Os[o_ind][1] == tuple(out_arr...)
                gyration += -f.memory + i % n
                break
            else
                circshift!(out_arr, 1)
            end
        end
    end

    return gyration % n
end

"alias of `gyration_number`"
gy(f, n) = gyration_number(f, n)

"""
    sign_gyration_compatibility_condition(f::BlockMap, n)

Compute the sign gyration compatibility condition of the given block map.
"""
function sign_gyration_compatibility_condition(f, n)
    to_return = gyration_number(f, n)
    k = n
    for i in 1:floor(Int, log2(n))
        if k % 2 != 0
            break
        else
            k = k ÷ 2
            to_return += os(f, k) * (n ÷ 2)
        end
    end

    return to_return % n
end

"alias of `sign_gyration_compatibility_condition`"
sgcc(f, n) = sign_gyration_compatibility_condition(f, n)


end # end of module BlockMaps