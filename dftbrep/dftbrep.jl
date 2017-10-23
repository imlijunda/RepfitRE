struct DftbRep

    nspl::Int
    cutoff::dp
    a1::dp
    a2::dp
    a3::dp

    knots::Vector{dp}
    coefs::Matrix{dp}

    function DftbRep(r::Vector{dp}, erep::Vector{dp}, s::dp = -1.0)

        knots, coefs = rep_coefs_3rd(r, erep; smoothness = s)
        r0 = knots[1]
        a1, a2, a3 = exp_coefs_3rd(coefs, r0)

        nspl = length(knots) - 1
        cutoff = knots[end]

        return new(nspl, cutoff, a1, a2, a3, knots, coefs)
    end

    function DftbRep(r_data::Vector{Vector{dp}}, e_data::Vector{Vector{dp}}, weight::Vector{dp}, smoothness::dp = -1.0)

        data = merge_scans(r_data, e_data, weight)
        if !check_scans(data)
            warn("Different repulsive energies were found at same distance.")
        end
        r, erep, w = resolve_scans(data)

        knots, coefs = rep_coefs_3rd(r, erep, w, smoothness)
        r0 = knots[1]
        a1, a2, a3 = exp_coefs_3rd(coefs, r0)

        nspl = length(knots) - 1
        cutoff = knots[end]


        return new(nspl, cutoff, a1, a2, a3, knots, coefs)
    end

    function DftbRep(skf::AbstractString)

        open(skf, "r") do f
            sdata = readlines(f)
        end

        f = x -> contains(x, "Spline")
        i = findfirst(f, sdata) + 1
        s = split(sdata[i])
        n = parse(Int, s[1])
        ptr_end = i + 1 + n

        data = sdata[i:ptr_end]

        s = split(data[1])
        nspl = parse(Int, s[1])
        cutoff = parse(dp, s[2])

        s = split(data[2])
        a1 = parse(dp, s[1])
        a2 = parse(dp, s[2])
        a3 = parse(dp, s[3])

        knots = zeros(dp, nspl + 1)
        coefs = zeros(dp, nspl, 6)
        for i = 1:nspl
            s = split(data[i+2])
            knots[i] = parse(dp, s[1])
            if i == nspl
                knots[i + 1] = parse(dp, s[2])
            end
            idx = 0
            for j = 3:length(s)
                idx += 1
                coefs[i, idx] = parse(dp, s[j])
            end
        end

        return new(nspl, cutoff, a1, a2, a3, knots, coefs)
    end
end

function exp_coefs_3rd(sp_coefs::Matrix{dp}, r0::dp)

    # Even today I still feel that the formalism of the exponential part is
    # dangerous and prone to domain error. Maybe we should explore other potential
    # form of the exponential part?
    c0 = sp_coefs[1, 1]
    c1 = sp_coefs[1, 2]
    c2 = sp_coefs[1, 3]
    #c3 is not needed
    a1 = -2.0 * c2 / c1
    a2 = log(-c1 / a1) + a1 * r0
    a3 = c0 - exp(-a1 * r0 + a2)

    return a1, a2, a3
end

function rep_coefs_3rd(r::Vector{dp}, erep::Vector{dp}; weight::Vector{dp} = Vector{dp}(0), smoothness::dp = -1.0)

    if isempty(weight)
        if smoothness < 0.0
            spl = Spline1D(r, erep; k=3, bc="extrapolate")
        else
            spl = Spline1D(r, erep; k=3, bc="extrapolate", s=smoothness)
        end
    else
        if smoothness < 0.0
            spl = Spline1D(r, erep; k=3, bc="extrapolate", w=weight)
        else
            spl = Spline1D(r, erep; k=3, bc="extrapolate", w=weight, s=smoothness)
        end
    end

    knots = get_knots(spl)
    #evaluate all values at knots
    n_knots = length(knots)
    coefs = zeros(dp, n_knots - 1, 6)
    #Spline value at ith knot
    Si = evaluate(spl, knots[1])
    Si_p = 0.0
    si_pp = 0.0
    #length of the spline
    R = 0.0

    for i = 1:n_knots - 1
        Si_p = derivative(spl, knots[i]; nu=1)
        Si_pp = derivative(spl, knots[i]; nu=2)
        coefs[i, 1] = Si        #0th order
        coefs[i, 2] = Si_p      #1st order
        coefs[i, 3] = Si_pp / 2 #2nd order
        #Solve 3rd order term at the end of current spline
        R = knots[i+1] - knots[i]
        Si = evaluate(spl, knots[i+1])
        coefs[i, 4] = (Si - coefs[i, 1] - coefs[i, 2] * R - coefs[i, 3] * R^2) / R^3
    end

    i = n_knots - 1
    #solve BC=zero linear system for last spline
    c0 = coefs[i, 1]
    c1 = coefs[i, 2]
    c2 = coefs[i, 3]
    R2 = R^2
    R3 = R^3
    R4 = R^4
    R5 = R^5
    LHS = reshape([  R3,    R4,    R5,
                   3*R2,  4*R3,  5*R4,
                   6*R , 12*R2, 20*R3 ], 3, 3).'
    RHS = [-(c0 + c1*R +   c2*R2),
           -(     c1   + 2*c2*R ),
           -(              c2   )]
    c345 = LHS\RHS;
    coefs[i, 4] = c345[1]
    coefs[i, 5] = c345[2]
    coefs[i, 6] = c345[3]

    return knots, coefs
end

default_weight(v::Vector{dp}) = ones(dp, length(v))
function merge_scans(r_data::Vector{Vector{dp}}, e_data::Vector{Vector{dp}}, weight::Vector{dp})

    total_size = 0
    for i = 1:length(r_data)
        total_size += length(r_data[i])
    end
    data = zeros(dp, total_size, 3)
    idx = 0
    for i = 1:length(r_data)
        for j = 1:length(r_data[i])
            idx += 1
            data[idx, 1] = r_data[i][j]
            data[idx, 2] = e_data[i][j]
            data[idx, 3] = weight[i]
        end
    end
    data = sortrows(data)

    return data
end

function check_scans(data::Matrix{dp}, threshold = PRECIS)

    (nrow, ncol) = size(data)
    for i = 2:nrow
        if equal(data[i-1, 1], data[i, 1]) && !equal(data[i-1, 2], data[i, 2], threshold)
            return false
        end
    end

    return true
end

function resolve_scans(data::Matrix{dp})

    r = dp[]
    erep = dp[]
    w = dp[]

    (nrow, ncol) = size(data)
    push!(r, data[1, 1])
    push!(erep, data[1, 2])
    push!(w, data[1, 3])
    for i = 2:nrow
        if !equal(r[end], data[i, 1])
            push!(r, data[i, 1])
            push!(erep, data[i, 2])
            push!(w, data[i, 3])
        end
    end

    return r, erep, w
end
