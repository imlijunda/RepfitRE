function vasp_esigma0(outcar::AbstractString)

    open(outcar, "r") do f
        data = readlines(f)
    end
    n = length(data)

    esigma0 = Inf
    for i = n:-1:1
        if startswith(strip(data[i]), "energy  without entropy=")
            s = split(data[i])
            esigma0 = parse(dp, s[7])
            break
        end
    end

    return esigma0
end
