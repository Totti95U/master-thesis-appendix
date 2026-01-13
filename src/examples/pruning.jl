using ProgressMeter

include("Pruning.jl")

hm = HenonMap(5.59, -1)
max_iter = 13

# TODO: [-5, 5] x [-5, 5] で多様体を区切って効率的に交点を探せるようにする
fixpt = hyperbolic_fixed_points(hm)[2]
smfds, umfds = manifolds(hm, num_iterations=max_iter, fixed_pt=fixpt, tolerance=1e-3)
println("Complete computing manifolds")

hmc_pts_info = intersections(smfds, umfds)
homoclinic_pts = Point2o[info[5] for info in hmc_pts_info]
println("Complete computing intersections")
println("Number of homoclinic points: ", length(homoclinic_pts))

### Symbolic encoding
threshold = 0.15
partition(pt::Point2o) = pt[1] < threshold ? 0 : 1
prim_sb, prim_ub = primary_branch(hm, partition, fixed_pt=fixpt)
symb_codes = symbolic_encoding(hm, homoclinic_pts, partition, prim_sb, prim_ub, iter=max_iter)
println("Complete symbolic encoding")

### Ploting field

# get default color palette
function plot_manifold()
    colors = Makie.current_default_theme()[:palette][:color] |> Makie.to_value

    fig = Figure(size=(500, 550))
    ax = Axis(fig[1, 1], aspect=1)
    vlines!(ax, [threshold], color=:black, linestyle=:dash)
    for smfd in smfds
        lines!(ax, smfd, color=colors[1])
        # scatter!(ax, smfd[1], color=colors[1], markersize=10)
    end
    for umfd in umfds
        lines!(ax, umfd, color=colors[2])
        # scatter!(ax, umfd[1], color=colors[2], markersize=10)
    end
    lines!(ax, prim_sb, color=:pink, linewidth=2, linestyle=:dash)
    lines!(ax, prim_ub, color=:green, linewidth=2, linestyle=:dash)
    scatter!(ax, homoclinic_pts, color=:red, markersize=8)

    xlims!(ax, -5.5, 5.5)
    ylims!(ax, -5.5, 5.5)

    return fig
end

function plot_symbolic_plane(bwd_div, fwd_div, symb_codes=symb_codes)
    fig = Figure(size=(500, 500))
    ax = Axis(fig[1, 1], aspect=1)

    fwd_div += 1

    println("need points at least ", 2^(bwd_div + fwd_div), " to fill the grid")
    if length(symb_codes) < 2^(bwd_div + fwd_div)
        @warn "Not enough homoclinic points to fill the grid"
    end

    xs = in_symbolic_plane.(symb_codes)
    for x in xs
        poly!(ax, [x, x, x, x] + [Point2o(0, 0), Point2o(1/2^fwd_div, 0), Point2o(1/2^fwd_div, 1/2^bwd_div), Point2o(0, 1/2^bwd_div)], color=:lightgray, strokewidth=0.5, strokecolor=:black)
    end
    # scatter!(ax, in_symbolic_plane.(symb_codes), color=:red, markersize=8)
    xlims!(ax, -0.05, 1.05)
    ylims!(ax, -0.05, 1.05)

    function lexi_to_knead(code::Vector{Int})
        w = similar(code)
        acc = accumulate(+, code)
        w[1] = code[1]
        for i in axes(w[1:end-1], 1)
            w[i+1] = sum(w[1:i]) % 2 == 0 ? code[i+1] : 1 - code[i+1]
        end
        return w
    end

    # generate all binary sequences of length `code_len`
    function keading_label(divnum::Int, flip=false)
        all_codes = Iterators.product(fill(0:1, divnum)...) |> collect
        all_codes = [collect(all_codes[i]) for i in 1:length(all_codes)] |> sort
        all_codes = lexi_to_knead.(all_codes)
        if flip
            all_codes = reverse.(all_codes)
        end

        # set axis ticks and labels
        ticks = (0:2^divnum-1) ./ 2^divnum .+ 1/(2^(divnum+1))
        tick_labels = [join(c) for c in all_codes]
        return ticks, tick_labels
    end
    
    ax.xticks = keading_label(fwd_div, false)
    ax.yticks = keading_label(bwd_div, true)

    # rotate x-tick labels by 90 degrees
    ax.xticklabelrotation = π/2

    return fig
end

function plot_primary_branches()
    fig1 = Figure(size=(400, 400))
    ax = Axis(fig1[1, 1], aspect=1)

    orbs = Point2o[]

    ub_seg, ub_tree = prepare_polyline(prim_ub)
    fwd_pts = copy(homoclinic_pts)
    fwd_d2s = similar(fwd_pts, MyFloat)
    @showprogress for i in eachindex(fwd_pts)
        best_d2 = typemax(MyFloat)
        pt = fwd_pts[i]
        for _ in 1:iter
            pt = forward(hm, pt)
            !isfinite(pt) && break
            d2, _, si, _ = closest_point(pt, ub_seg, ub_tree)
            if d2 < best_d2
                best_d2 = d2
                fwd_pts[i] = pt
            end
            push!(orbs, pt)
        end
        fwd_d2s[i] = best_d2
    end

    sb_seg, sb_tree = prepare_polyline(prim_sb)
    bwd_pts = copy(homoclinic_pts)
    bwd_d2s = similar(bwd_pts, MyFloat)
    @showprogress for i in eachindex(bwd_pts)
        best_d2 = typemax(MyFloat)
        pt = bwd_pts[i]
        for _ in 1:iter
            pt = backward(hm, pt)
            !isfinite(pt) && break
            d2, _, si, _ = closest_point(pt, sb_seg, sb_tree)
            if d2 < best_d2
                best_d2 = d2
                bwd_pts[i] = pt
            end
            push!(orbs, pt)
        end
        bwd_d2s[i] = best_d2
    end

    lines!(ax, prim_sb, color=:blue, linewidth=2, linestyle=:dash)
    lines!(ax, prim_ub, color=:red, linewidth=2, linestyle=:dash)

    scatter!(ax, homoclinic_pts, color=:gray, markersize=6)

    scatter!(ax, fwd_pts, color=:blue, markersize=8, marker=:xcross)
    scatter!(ax, bwd_pts, color=:red, markersize=8, marker=:xcross)

    println("min fwd d2: ", minimum(fwd_d2s))
    println("max fwd d2: ", maximum(fwd_d2s))
    println("min bwd d2: ", minimum(bwd_d2s))
    println("max bwd d2: ", maximum(bwd_d2s))

    return fig1
end

function plot_primary_pruned_region(p_blocks, bwd_div=3, fwd_div=3)
    fig = Figure(size=(400, 400))
    ax = Axis(fig[1, 1], aspect=1)
    xlims!(ax, -0.05, 1.05)
    ylims!(ax, -0.05, 1.05)
    
    for block in p_blocks
        a, b, c, d = in_symbolic_plane.(block)
        poly!(ax, [a, b, c, d], color=:red, strokewidth=0.5, strokecolor=:black)
    end

    scatter!(ax, in_symbolic_plane.(symb_codes), color=:blue, markersize=4)

    function lexi_to_knead(code::Vector{Int})
        w = similar(code)
        acc = accumulate(+, code)
        w[1] = code[1]
        for i in axes(w[1:end-1], 1)
            w[i+1] = sum(w[1:i]) % 2 == 0 ? code[i+1] : 1 - code[i+1]
        end
        return w
    end

    function keading_label(divnum::Int, flip=false)
        all_codes = Iterators.product(fill(0:1, divnum)...) |> collect
        all_codes = [collect(all_codes[i]) for i in 1:length(all_codes)] |> sort
        all_codes = lexi_to_knead.(all_codes)
        if flip
            all_codes = reverse.(all_codes)
        end

        # set axis ticks and labels
        ticks = [in_symbolic_plane(HomoclinicCode(Int[], c))[1] for c in all_codes]
        tick_labels = [join(c) for c in all_codes]
        return ticks, tick_labels
    end
    
    ax.xticks = keading_label(fwd_div, false)
    ax.yticks = keading_label(bwd_div, true)
    ax.xticklabelrotation = π/2

    return fig
end

blocks = primary_pruning_front(symb_codes)
# TODO: block から forbidden word を抽出する