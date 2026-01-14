using Graphs

const Word = Vector{Int}

_word_to_int(w::Word) = sum(w .* (2 .^ (0:length(w)-1))) + 1
function _int_to_word(x::Int; pad::Int=1)
    w = Word()
    for _d in reverse(string(x, base=2, pad=pad))
        d = parse(Int, _d)
        push!(w, d)
    end
    return w
end

"""
    in_shift(w::Word, B::Vector{Word})

Check if word `w` lies in the shift defined by forbidden words `B`.
"""
function in_shift(w::Word, B::Vector{Word})
    # full shift case
    if isempty(B)
        return true
    end

    m = maximum(length.(B))

    # generate transition matrix
    A = zeros(Int, (2^(m-1), 2^(m-1)))
    for i in axes(A, 1)
        A[i, i] = 1
    end
    # A is indexed by ordering like 0000, 1000, 0100, 1100, 0010, 1010, ..., 0001, 1001, 0101, 1101, 0011, 1011, ..., 1111
    A = repeat(A, inner=(2, 1), outer=(1, 2))

    # fill row and column for each forbidden word with 0s
    ind_to_del = [_word_to_int(b) for b in B]
    for ind in ind_to_del
        A[ind, :] .= 0
        A[:, ind] .= 0
    end

    # find words in strong connected components with cycles
    cyclic_comps = Vector{Vector{Int}}()
    cyclic_vertices = Set{Int}()
    gA = DiGraph(A)
    scc = strongly_connected_components(gA)
    for comp in scc
        if length(comp) > 1 || A[comp[1], comp[1]] == 1
            push!(cyclic_comps, comp)
            union!(cyclic_vertices, comp)
        end
    end

    # test w
    if length(w) <= m
        # check if w lies in SCC
        return (_word_to_int(w) in cyclic_vertices) ? true : false
    else
        # check every length-m subword of w lies in V
        for i in 1:length(w)-m+1
            if _word_to_int(w[i:i+m-1]) in ind_to_del
                return false
            end
        end

        # check w can come from cyclic_words
        can_comes_from = false
        for comp in cyclic_comps
            if has_path(gA, comp[1], _word_to_int(w[1:m]))
                can_comes_from = true
                break
            end
        end
        if !can_comes_from
            false
        end

        # check w can go to cyclic_words
        can_goes_to = false
        for comp in cyclic_comps
            if has_path(gA, _word_to_int(w[end-m+1:end]), comp[1])
                can_goes_to = true
                break
            end
        end
        return can_goes_to
    end
end

"""
    B1 ⊏ B2

Check if the shift defined by forbidden words `B1` is contained in the shift defined by forbidden words `B2`.
"""
function ⊏(B1::Vector{Word}, B2::Vector{Word})
    for b in B2
        in_shift(b, B1) && return false
    end
    return true
end

"""
    B1 ⊐ B2

Check if the shift defined by forbidden words `B1` contain the shift defined by forbidden words `B2`.
"""
function ⊐(B1::Vector{Word}, B2::Vector{Word})
    return B2 ⊏ B1
end

"""
    shift_hasse(V::Vector{Vector{Word}})

Construct the Hasse diagram of the subshift lattice defined by the list of forbidden word sets `V`.
"""
function shift_hasse_edge(V::Vector{Vector{Word}})
    n = length(V)
    # construct all edges
    E = Vector{Tuple{Int,Int}}()
    for i in axes(V, 1), j in axes(V, 1)
        if i != j && V[i] ⊏ V[j]
            push!(E, (i, j))
        end
    end

    # remove transitive edges
    keep = trues(length(E))
    for (k, (u, v)) in enumerate(E)
        for l in axes(V, 1)
            if (u, l) in E && (l, v) in E
                keep[k] = false
                break
            end
        end
    end
    return E[keep]
end

"""
    shift_hasse_diagram(V)

与えられた SFT の間の包含関係によるハッセ図を出力します.

`V` の各元は SFT の forbidden word のリストです.
図内の頂点の番号は V のインデックスに対応しています
"""
function shift_hasse_diagram(V::Vector{Vector{Word}})
    edges = shift_hasse_edge(V)
    n = length(V)
    G = SimpleDiGraph(n)
    for (i, j) in edges
        add_edge!(G, i, j)
    end

    # compute rank (levels)
    rank = fill(0, n)
    for v in topological_sort(G)
        preds = inneighbors(G, v)
        if !isempty(preds)
            rank[v] = 1 + maximum(rank[p] for p in preds)
        end
    end

    # layout positions by rank
    levels = Dict(r => findall(x -> rank[x] == r, 1:n) for r in unique(rank))
    positions = zeros(2, n)
    for (r, nodes) in sort(collect(levels))
        added_node = Int[]
        for (k, v) in enumerate(nodes)
            tmp_x = isempty(inneighbors(G, v)) ? 0.0 : minimum(positions[1, inneighbors(G, v)])
            for u in added_node
                if any(abs.(positions[1, added_node] .- tmp_x) .< 1.0)
                    tmp_x += 1.0
                else
                    break
                end
            end
            push!(added_node, v)
            positions[1, v] = tmp_x
            positions[2, v] = -r
        end
    end
    # 分岐ノードを真ん中にする
    # for (r, nodes) in sort(collect(levels))
    #     for v in nodes
    #         if length(outneighbors(G, v)) > 1
    #             child_xs = positions[1, outneighbors(G, v)] |> unique
    #             positions[1, v] = (maximum(child_xs) + minimum(child_xs)) / 2
    #         end
    #     end
    # end
    # 分岐を受けるノードも真ん中にする
    # for (r, nodes) in sort(collect(levels); rev=true)
    #     for v in nodes
    #         if length(inneighbors(G, v)) > 1
    #             parent_xs = positions[1, inneighbors(G, v)] |> unique
    #             positions[1, v] = (maximum(parent_xs) + minimum(parent_xs)) / 2
    #         end
    #     end
    # end

    # draw
    gplot(G, positions[1, :], positions[2, :]; nodelabel=1:n)
end
