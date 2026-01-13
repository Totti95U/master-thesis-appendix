include("elem-arith.jl")

# counting the number of periodic points with primitive method
function count_fix(markers::Vector, n::Int)
    # initialize the counter
    counter = 0

    # generate all words constructed with 1 and 2
    for i in 1:2^n
        # decode `i` into the binary sequence `b_seq`
        b_seq = Int[]
        for _ in 1:n
            push!(b_seq, i % 2 + 1)
            i = i >> 1
        end

        for marker in markers, j in 1:n
            l = length(marker)

            # split the marker into two vector
            marker1 = [s == "*" ? 1 : s for s in marker]
            marker2 = [s == "*" ? 2 : s for s in marker]

            # extend `b_seq` for looping boundary
            # the length must have
            must_l = l + n - 1
            ext_b_seq = repeat(b_seq, must_l÷n)
            ext_b_seq = [ext_b_seq; b_seq[1:must_l%n]]
            
            # identify the `b_seq` has `marker` as subvector
            if ext_b_seq[j:j+l-1] == marker1 || ext_b_seq[j:j+l-1] == marker2
                counter += 1
                break
            end
        end
    end

    return counter
end
"alias for count_fix"
count_fix(marker::Tuple, n::Int) = count_fix([marker], n)

# count the number of (2, n)-cascade for all n ≤ m
function count_cascade(markers::Vector, m::Int)
    fn = Int[]
    pn = Int[]
    on = Int[]
    cascade_n = Int[]
    for n in 1:m
        push!(fn, count_fix(markers, n))
        push!(pn, sum(d -> mobius(d) * fn[Int(n/d)], divisors(n)))
        push!(on, pn[end]÷n)
        b_d = binary_divisors(n)
        push!(cascade_n, on[end])
        if length(b_d) > 0
            cascade_n[end] -= sum([cascade_n[d] for d in b_d])
        end
        cascade_n[end] ÷= 2
        # print(cascade_n[end],", ")
    end
    # println()
    return cascade_n
end
"alias for count_cascade"
count_cascade(marker::Tuple, m::Int) = count_cascade([marker], m)

function gcd2_cascade(markers::Vector, m::Int)
    fn = Int[]
    cascade_n = Int[]
    for n in 1:m
        push!(fn, count_fix(markers, n))
        c_n = 0
        for d in divisors(n)
            c_n += gcd(2, d) * mobius(d) * fn[Int(n/d)]
        end
        push!(cascade_n, c_n ÷ 2n)
    end
    return cascade_n
end
"alias for gcd2_cascade"
gcd2_cascade(marker::Tuple, m::Int) = gcd2_cascade([marker], m)

"the number of periodic orbit with even1's"
function _even1_orbit(markers::Vector, n::Int)
    # initialize the counter
    count = 0

    # generate all words constructed with 1 and 2
    for i in 1:2^n
        # decode `i` into the binary sequence `b_seq`
        b_seq = Int[]
        for _ in 1:n
            push!(b_seq, i % 2 + 1)
            i = i >> 1
        end

        p = period(b_seq)
        if p < n
            continue
        end

        for marker in markers, j in 1:n
            l = length(marker)

            # split the marker into two vector
            marker1 = [s == "*" ? 1 : s for s in marker]
            marker2 = [s == "*" ? 2 : s for s in marker]

            # extend `b_seq` for looping boundary
            # the length must have
            must_l = l + n - 1
            ext_b_seq = repeat(b_seq, must_l÷n)
            ext_b_seq = [ext_b_seq; b_seq[1:must_l%n]]
            
            # identify the `b_seq` has `marker` as subvector
            if ext_b_seq[j:j+l-1] == marker1 || ext_b_seq[j:j+l-1] == marker2
                count += (sum(b_seq) % 2 == 0)
                break
            end
        end
    end

    return count ÷ n
end

even1_cascade(markers::Vector, m::Int) = map(n -> _even1_orbit(markers, n), 1:m)
even1_cascade(marker::Tuple, m::Int) = even1_cascade([marker], m)
