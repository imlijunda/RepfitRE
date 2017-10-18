function eval_rep(rep::DftbRep, r::dp)

    erep = 0.0
    if r < rep.knots[1]
        erep = exp(-rep.a1 * r + rep.a2) + rep.a3
        return erep
    elseif r > rep.cutoff
        erep = 0.0
        return erep
    else
        for i = rep.nspl:-1:1
            if r >= rep.knots[i]
                coefs = rep.coefs[i, :]
                r0 = rep.knots[i]
                for j = 1:length(coefs)
                    erep += coefs[j] * ((r - r0)^(j - 1))
                end
                return erep
            end
        end
    end

    return Inf
end

function eval_rep(rep::DftbRep, r::Vector{dp})

    n = length(r)
    erep = zeros(dp, n)
    for i = 1:n
        erep[i] = eval_rep(rep, r[i])
    end

    return erep
end

function rms_rep(rep1::DftbRep, rep2::DftbRep; start=0.01, np=1000)

    cutoff = min(rep1.cutoff, rep2.cutoff)

    s = 0.0
    r = linspace(start, cutoff, np)
    for i = 1:np
        s = (eval_rep(rep1, r[i]) - eval_rep(rep2, r[i]))^2.0
    end
    s = sqrt(s / np)

    return s
end

function smooth_rep(rep::DftbRep, n::Int, smoothness::dp)

    r0 = rep.knots[1]
    cutoff = rep.cutoff

    r = linspace(r0, cutoff, n)
    r = collect(r)
    erep = eval_rep(rep, r)

    new_rep = DftbRep(r, erep; s = smoothness)

    return new_rep
end

function write_rep(filepath::AbstractString, rep::DftbRep)

    open(filepath, "w") do f
        b = 0
        b += write(f, "Spline\n")

        s = @sprintf("%d  %.15f\n", rep.nspl, rep.cutoff)
        b += write(f, s)

        s = @sprintf("%.15f  %.15f  %.15f\n", rep.a1, rep.a2, rep.a3)
        b += write(f, s)

        for i = 1:rep.nspl - 1
            s = @sprintf("%.15f  %.15f  %.15f  %.15f  %.15f  %.15f\n",
                    rep.knots[i], rep.knots[i+1], rep.coefs[i, 1:4]...)
            b += write(f, s)
        end
        i = rep.nspl
        s = @sprintf("%.15f  %.15f  %.15f  %.15f  %.15f  %.15f  %.15f  %.15f\n",
                rep.knots[i], rep.knots[i+1], rep.coefs[i, :]...)
        b += write(f, s)
    end

    return b
end
