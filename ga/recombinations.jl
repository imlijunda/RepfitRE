function discrete(c1::Vector{T}, c2::Vector{T}) where {T}

    n = length(c1)
    assert(n == length(c2))
    c  = zeros(T, n)

    for i = 1:n
        c[i] = prob(0.5) ? c1[i]:c2[i]
    end

    return c
end

function intermediate(c1::Vector{T}, c2::Vector{T}, δ::T = 0.25) where {T<:Real}

    n1, n2 = length(c1), length(c2)
    assert(n1 == n2)

    #According to Julia doc, random numbers are already generated uniformly in [0, 1)
    α = rand(T, n1)
    #map to [-δ, 1 + δ]
    k::T = 1.0 + δ * 2
    b::T = -δ
    for i = 1:n1
        α[i] = α[i] * k + b
    end

    c = zeros(T, n1)
    for i = 1:n1
        c[i] = c1[i] + α[i] * (c2[i] - c1[i])
    end

    return c
end

function dynamic_intermediate(c1::Vector{T}, c2::Vector{T}, rm::Vector{T}, rp::Vector{T}) where {T<:Real}

    n1, n2, nr1, nr2 = length(c1), length(c2), length(rm), length(rp)
    assert(n1 == n2 == nr1 == nr2)

    α = rand(T, n1)
    for i = 1:n1
        r = rp[i] - rm[i]
        δm = (c2[i] - rm[i]) / r
        δp = (rp[i] - c1[i]) / r
        k::T = 1.0 + δm + δp
        b = -δm
        α[i] = a[i] * k + b
    end

    c = zeros(T, n1)
    for i = 1:n1
        c[i] = c1[i] + α[i] * (c2[i] - c1[i])
    end

    return c
end

function line(c1::Vector{T}, c2::Vector{T}, δ::T = 0.25) where {T<:Real}

    n1, n2 = length(c1), length(c2)
    assert(n1 == n2)

    #same as intermediate but only apply single α value
    α = rand(T)
    k::T = 1.0 + δ * 2
    b::T = -δ
    α = α * k + b

    c = zeros(T, n1)
    for i = 1:n1
        c[i] = c1[i] + α * (c2[i] - c1[i])
    end

    return c
end
