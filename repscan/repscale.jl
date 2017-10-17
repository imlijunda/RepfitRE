mutable struct RepScale

    rep::Repulsive

    d::dp
    w::dp
    # energies
    edft::dp
    eband::dp
    # total repulsive energy
    tot_erep::dp
    # reduced repulsive energy
    net_erep::dp

    scaled::Vector{RepScale}

    function RepScale(rep::Repulsive, sfun::Function)

        rsc = new()
        rsc.rep = rep

        rsc.d = rep.distance[1]
        rsc.w = rep.weight[1]
        rsc.edft = 0.0
        rsc.eband = 0.0
        rsc.tot_erep = 0.0
        rsc.net_erep = 0.0
        rsc.scaled = RepScale[]

        # scale repulsive if needed
        for i = 2:length(rep.distance)
            scale_ratio = rep.distance[i] / rsc.d
            newgeo = sfun(rep.nei.geo, scale_ratio)
            newnei = Neighbour(newgeo, rep.cutoff)
            newrep = Repulsive(newnei, rep.sp1, rep.sp2, rep.precision)
            push!(rsc.scaled, RepScale(newrep, sfun))
        end

        return rsc
    end
end

function uni_scale(geo::Geometry, r::dp)

    newgeo = Geometry(geo)
    newgeo.lattice = newgeo.lattice .* r

    return newgeo
end

function uni_max_gen(rep::Repulsive)

    if length(rep.distance) < 2
        return 1
    end

    k = 0
    d = rep.distance[1]
    r = rep.distance[2] / d
    while true
        k += 1
        d *= r
        if d >= rep.cutoff
            break
        end
    end

    return k
end

function uni_max_cutoff(rep::Repulsive, gen::Int)

    if length(rep.distance) < 2
        return 0.0
    end

    d = rep.distance[1]
    r = rep.distance[2] / d
    for i = 1:gen
        d *= r
    end

    return d
end

function shift_tot_erep!(rsc::RepScale, ee::dp)

    rsc.tot_erep -= ee
    for i = 1:length(rsc.scaled)
        shift_tot_erep!(rsc.scaled[i], ee)
    end

    return rsc
end

function calc_net_erep!(rsc::RepScale)

    tot = 0.0
    for i = 1:length(rsc.scaled)
        calc_net_erep!(rsc.scaled[i])
        tot += rsc.scaled[i].net_erep * rsc.scaled[i].w
    end
    rsc.net_erep = (rsc.tot_erep - tot) / rsc.w

    return rsc
end
