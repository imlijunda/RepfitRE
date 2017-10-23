function truncation(fitness::Vector, n::Int)

    p = sortperm(fitness; rev = true)

    return p[1:n]
end

function __select(p::Vector, n::Int)

    selected = Vector{Int}(n)
    s = cumsum(p)
    for i = 1:n
        r = rand()
        idx = findfirst(x -> x >= r, s)
        selected[i] = idx
    end

    return selected
end

function roulette(fitness::Vector, n::Int)

    p = fitness ./ sum(fitness)

    return __select(p, n)
end


function __ranking_linear(fitness::Vector, n::Int, ηm::dp, ηp::dp)

    N = length(fitness)
    p = Vector{dp}(N)
    ranks = sortperm(fitness)
    for i = 1:N
        p[i] = (ηm + (ηp - ηm) * (ranks[i] - 1.0) / (N - 1.0)) / N
    end

    return __select(p, n)
end
function ranking_linear(fitness::Vector, n::Int, η::dp)

    assert(0.0 <= η <= 2.0)
    ηm = η
    ηp = 2.0 - η

    return __ranking_linear(fitness, n, ηm, ηp)
end
function ranking_linear_gradient(fitness::Vector, n::Int, rm::dp)

    ηm = 2.0 / (rm + 1.0)
    ηp = 2.0 * rm / (rm + 1.0)

    return __ranking_linear(fitness, n, ηm, ηp)
end

function ranking_exp(fitness::Vector, n::Int, c::dp)

    assert(0.0 < c < 1.0)

    N = length(fitness)
    p = Vector{dp}(N)
    ranks = sortperm(fitness)

    coefs = (c - 1.0) / (c^N - 1.0)
    for i = 1:N
        p[i] = coefs * c^(N - ranks[i])
    end

    return __select(p, n)
end

function tournament(fitness::Vector, n::Int, gp_size::Int)

    ranks = sortperm(fitness)
    N = length(fitness)
    winners = Int[]
    for i = 1:n
        candidates = randperm(N)[1:gp_size]
        f = fitness[candidates]
        m = 0.0
        best = 0
        for j = 1:gp_size
            if f[j] > m
                best = candidates[j]
                m = f[j]
            end
        end
        push!(winners, best)
    end

    return winners
end
