scan_vector(start::Vector{dp}, stop::Vector{dp}) = stop - start

mutable struct LinearScan

    vect::Vector{dp}
    n::Int
    delta::dp

    spfunc::Function

    function LinearScan(v::Vector{dp}, n::Int; spfunc::Function = linspace)

        ls = new()
        ls.vect = copy(v)
        ls.n = n
        ls.delta = 0.0

        ls.spfunc = spfunc

        return ls
    end
end

function shift_scan!(ls::LinearScan, distance::dp)

    ls.delta = distance

    return ls
end

function coord_list(ls::LinearScan)

    r = ls.spfunc(0.0, 1.0, ls.n)
    r += ls.delta

    clist = Vector{Vector{dp}}(ls.n)
    for i = 1:ls.n
        clist[i] = ls.vect * r[i]
    end

    return clist
end

function geo_list(ls::LinearScan, geo::Geometry, spart::AtomList)

    clist = coord_list(ls)
    glist = Vector{Geometry}(ls.n)

    for i = 1:ls.n
        newgeo = Geometry(geo)
        move!(newgeo, spart, clist[i])
        glist[i] = newgeo
    end

    return glist
end
