mutable struct Neighbour

    #A pointer to refered geometry
    geo::Geometry

    cutoff::dp
    #Number of neighbours for ith atom
    nnei:: Vector{Int}
    #Index of ith atom's jth neighbour
    inei::Vector{Vector{Int}}
    #Distances
    dnei::Vector{Vector{dp}}

    function Neighbour(geo::Geometry, cutoff::dp)

        nei = new()
        nei.geo = geo
        update_cutoff!(nei, cutoff)

        return nei
    end

end

function update_cutoff!(nei::Neighbour, cutoff::dp)

    geo = nei.geo
    natom = size(geo.coord, 2)

    nnei = zeros(Int, natom)
    inei = Vector{Vector{Int}}(natom)
    dnei = Vector{Vector{dp}}(natom)
    for i = 1:natom
        inei[i] = Vector{Int}()
        dnei[i] = Vector{dp}()
    end

    # search ranges of shifted lattices
    lpts = lat_pts(geo.lattice, cutoff)
    nlpts = size(lpts, 2)
    #distance map
    m = zeros(dp, natom, natom)
    for l = 1:nlpts
        dmap!(m, geo.lattice, geo.coord, @view(lpts[:, l]))
        #Ref: DFTB+ code
        for i = 1:natom
            for j = 1:i
                if (m[i, j] <= cutoff) && (m[i, j] > 0)
                    nnei[i] += 1
                    push!(inei[i], j)
                    push!(dnei[i], m[i, j])
                end
            end
        end
    end

    nei.cutoff = cutoff
    nei.nnei = nnei
    nei.inei = inei
    nei.dnei = dnei

    return nei
end

# Ref: DFTB+ code
function lat_pts(lattice::Matrix{dp}, cutoff::dp, neg_ext::Int = 1, pos_ext::Int = 1)

    iv = reciprocal(lattice)

    # search ranges
    tmp = floor(Int, cutoff * norm(iv[:, 1]))
    rngA = -tmp - neg_ext : tmp + pos_ext
    szA = length(rngA)
    tmp = floor(Int, cutoff * norm(iv[:, 2]))
    rngB = -tmp - neg_ext : tmp + pos_ext
    szB = length(rngB)
    tmp = floor(Int, cutoff * norm(iv[:, 3]))
    rngC = -tmp - neg_ext : tmp + pos_ext
    szC = length(rngC)
    sz = szA * szB * szC

    lat_points = zeros(Int, 3, sz);
    idx = 0;
    for i = rngA
        for j = rngB
            for k = rngC
                idx += 1
                lat_points[1, idx] = i
                lat_points[2, idx] = j
                lat_points[3, idx] = k
            end
        end
    end

    return lat_points
end

function dmap(lattice::Matrix{dp}, coord::Matrix{dp}, lat_shift::Union{Vector{Int}, SubArray{Int, 1}})

    n = size(coord, 2)
    abs_shift = lattice * lat_shift
    abs_coord = lattice * coord

    m = zeros(dp, n, n)
    delta = zeros(dp, 3)
    for i = 1:n
        for j = 1:n
            delta[1] = abs_coord[1, i] - (abs_coord[1, j] + abs_shift[1])
            delta[2] = abs_coord[2, i] - (abs_coord[2, j] + abs_shift[2])
            delta[3] = abs_coord[3, i] - (abs_coord[3, j] + abs_shift[3])
            m[i, j] = norm(delta)
        end
    end

    return m
end
function dmap!(m::Matrix{dp}, lattice::Matrix{dp}, coord::Matrix{dp}, lat_shift::Union{Vector{Int}, SubArray{Int, 1}})

    n = size(coord, 2)
    abs_shift = lattice * lat_shift
    abs_coord = lattice * coord

    delta = zeros(dp, 3)
    for i = 1:n
        for j = 1:n
            delta[1] = abs_coord[1, i] - (abs_coord[1, j] + abs_shift[1])
            delta[2] = abs_coord[2, i] - (abs_coord[2, j] + abs_shift[2])
            delta[3] = abs_coord[3, i] - (abs_coord[3, j] + abs_shift[3])
            m[i, j] = norm(delta)
        end
    end

    return m
end

function dmap(coord1::Matrix{dp}, coord2::Matrix{dp})

    (s1, n1) = size(coord1)
    (s2, n2) = size(coord2)
    @assert s1 == s2 == 3

    delta = zeros(dp, 3)
    m = zeros(dp, n1, n2)
    for i = 1:n1
        for j = 1:n2
            delta[1] = coord1[1, i] - coord2[1, j]
            delta[2] = coord1[2, i] - coord2[2, j]
            delta[3] = coord1[3, i] - coord2[3, j]
            m[i, j] = norm(delta)
        end
    end

    return m
end
function dmap!(m::Matrix{dp}, coord1::Matrix{dp}, coord2::Matrix{dp})

    (s1, n1) = size(coord1)
    (s2, n2) = size(coord2)
    @assert s1 == s2 == 3
    @assert size(m) == (n1, n2)

    delta = zeros(dp, 3)
    for i = 1:n1
        for j = 1:n2
            delta[1] = coord1[1, i] - coord2[1, j]
            delta[2] = coord1[2, i] - coord2[2, j]
            delta[3] = coord1[3, i] - coord2[3, j]
            m[i, j] = norm(delta)
        end
    end

    return m
end
