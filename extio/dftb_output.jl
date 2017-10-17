function dftb_mermin(detailed_out::AbstractString)

    open(detailed_out, "r") do f
        data = readlines(f)
    end
    n = length(data)

    mermin = Inf
    for i = n:-1:1
        if startswith(strip(data[i]), "Total Mermin free energy")
            s = split(data[i])
            mermin = parse(dp, s[5])
            break
        end
    end

    return mermin
end

function dftb_repulsive(detailed_out::AbstractString)

    open(detailed_out, "r") do f
        data = readlines(f)
    end
    n = length(data)

    repulsive = Inf
    for i = n:-1:1
        if startswith(strip(data[i]), "Repulsive energy")
            s = split(data[i])
            repulsive = parse(dp, s[3])
            break
        end
    end

    return repulsive
end
