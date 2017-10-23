# range: 0.1~0.5 * search invterval
function discrete_BGA!(chromosome::Vector{T}, range::Vector{T}; k = 16) where {T<:Real}

    #Calculate δ
    δ = zeros(T, k)
    pδ = 1.0 / k
    for i = 1:k
        δ[i] = prob(pδ) ? 0.0:2.0^(-i + 1)
    end
    for i = 1:length(chromosome)
        chromosome[i] += prob(0.5) ? -δ[i] * range[i] : δ[i] * range[i]
    end

    return chromosome
end

function continuous_BGA!(chromosome::Vector{T}, range::Vector{T}; k = 16) where {T<:Real}

    for i = 1:length(chromosome)
        δ::T = 2.0^(-k * rand(T))
        chromosome[i] += prob(0.5) ? -δ * range[i] : δ * range[i]
    end

    return chromosome
end

function permutation!(chromosome::Vector)

    n = length(chromosome)
    perm = randperm(n)
    chromosome = chromosome[perm]

    return chromosome
end
# variants of permutation
function shiftleft!(chromosome::Vector, n::Int = 0)

    len = length(chromosome)
    if n == 0
        n = rand(1:len)
    end

    idx = n
    tmp = copy(chromosome)
    for i = 1:len
        idx = rem(idx, len) + 1
        chromosome[i] = tmp[idx]
    end

    return chromosome
end
function shiftright!(chromosome::Vector, n::Int = 0)

    len = length(chromosome)
    if n == 0
        n = rand(1:len)
    end

    idx = n
    tmp = copy(chromosome)
    for i = 1:len
        idx = rem(idx, len) + 1
        chromosome[idx] = tmp[i]
    end

    return chromosome
end

# simply flip random bit(s)
function flip!(chromosome::Vector{Bool}; count = 1)

    n = length(chromosome)
    for i = 1:count
        ptr = rand(1:n)
        chromosome[ptr] = !chromosome[ptr]
    end

    return chromosome
end
