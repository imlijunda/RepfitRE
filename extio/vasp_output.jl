function vasp_esigma0(outcar::AbstractString)

    f=open(outcar, "r")
    data = readlines(f)
    n = length(data)
    close(f)

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

function vasp_esigma0_ha(outcar::AbstractString)

    ev = vasp_esigma0(outcar)

    return  ev * hartree_per_ev
end
