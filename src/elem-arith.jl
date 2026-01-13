# mobius function
function mobius(n::Int)
    if n < 1
        return 0
    end
    if n == 1
        return 1
    end
    # number of prime factors
    p = 0
    for i in 2:n
        if n % i == 0
            p += 1
            n ÷= i
        end
        if n % i == 0
            return 0
        end
    end
    return p % 2 == 0 ? 1 : -1
end

# divisor iterator
function divisors(n::Int)
    divs = Int[]
    for i in 1:floor(Int, sqrt(n))
        if n % i == 0
            push!(divs, i)
            if i != n ÷ i
                push!(divs, n ÷ i)
            end
        end
    end
    return sort!(divs)
end

# divide by power of 2
function binary_divisors(n::Int)
    divs = Int[]
    while n > 1
        if n % 2 != 0
            break
        end
        n ÷= 2
        push!(divs, n)
    end
    return divs
end

# return the period
function period(seq::Vector{Int})
    n = length(seq)
    for d in divisors(n)
        if all(seq[i] == seq[i + d] for i in 1:(n - d))
            return d
        end
    end
    return n
end