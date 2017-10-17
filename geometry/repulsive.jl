mutable struct Repulsive

    nei::Neighbour

    cutoff::dp
    sp1::Int
    sp2::Int

    precision::dp
    distance::Vector{dp}
    weight::Vector{dp}

    nrep::Int
    raw_rep::Matrix{dp}

    function Repulsive(nei::Neighbour, sp1::Int, sp2::Int, precision::dp = PRECIS)

        cutoff = nei.cutoff
        rep = new()
        rep.nei = nei
        rep.cutoff = cutoff
        update_raw!(rep, sp1, sp2)
        update_pdw!(rep, precision)

        return rep
    end

end

function update_cutoff!(rep::Repulsive, cutoff::dp)

    update_cutoff!(nei, cutoff)
    rep.cutoff = cutoff
    update_raw!(rep, rep.sp1, rep.sp2)
    update_pdw!(rep, rep.precision)

    return rep
end

function update_raw!(rep::Repulsive, sp1::Int, sp2::Int)

    sp = rep.nei.geo.sp
    natom = size(rep.nei.geo.coord, 2)
    ntmp = sum(rep.nei.nnei)
    tmp = zeros(dp, ntmp, 2)
    nrep = 0
    for i = 1:natom
        for j = 1:rep.nei.nnei[i]
            if ((sp[i] == sp1) && (sp[rep.nei.inei[i][j]] == sp2)) ||
                ((sp[i] == sp2) && (sp[rep.nei.inei[i][j]] == sp1))
                nrep += 1
                tmp[nrep, 1] = rep.nei.dnei[i][j]
                #in case of self-image
                tmp[nrep, 2] = rep.nei.inei[i][j] == i ? 0.5:1.0
            end
        end
    end
    raw_rep = sortrows(tmp[1:nrep, :])

    rep.sp1 = sp1
    rep.sp2 = sp2
    rep.nrep = nrep
    rep.raw_rep = raw_rep

    return rep
end

function update_pdw!(rep::Repulsive, precision::dp)

    distance = dp[]
    weight = dp[]
    if rep.nrep >0
        push!(distance, rep.raw_rep[1, 1])
        push!(weight, rep.raw_rep[1, 2])
        for i = 2:rep.nrep
            if !equal(distance[end], rep.raw_rep[i, 1], precision)
                push!(distance, rep.raw_rep[i, 1])
                push!(weight, rep.raw_rep[i, 2])
            else
                weight[end] += rep.raw_rep[i, 2]
            end
        end
    end

    rep.precision = precision
    rep.distance = distance
    rep.weight = weight

    return rep
end
