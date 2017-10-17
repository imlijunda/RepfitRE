#Ref: Simon K. Kearsley, Acta Cryst. (1989). A45, 208-210
function quaternion!(coord1::Matrix{dp}, coord2::Matrix{dp}, ipair::Matrix{Int}, weight::Vector{dp}, minimize = true)

    #recenter two molecules
    coord_center!(coord1, coord2, ipair, weight)
    #rotate second molecules for maximum overlap
    m = symm_matrix(coord1, coord2, ipair, weight)
    λ, Q = eig(m)
    # λ and Q are already sorted ascending
    # 1st row for maximized displacement, 4th for minimum
    row = minimize ? 4 : 1
    q = Q[row, :]
    rot = rot_matrix(q)
    new_coord = rot * coord2
    (n1, n2) = size(coord2)
    for i = 1:n1
        for j = 1:n2
            coord2[i, j] = new_coord[i, j]
        end
    end

    return coord1, coord2
end

function coord_rms(coord1::Matrix{dp}, coord2::Matrix{dp}, ipair::Matrix{Int}, weight::Vector{dp})

    n = length(weight)
    dx = 0.0
    dy = 0.0
    dz = 0.0
    dist2 = 0.0
    total = 0.0
    rms = 0.0
    for i = 1:n
        at1 = ipair[1, i]
        at2 = ipair[2, i]
        dx = coord1[1, i] - coord2[1, i]
        dy = coord1[2, i] - coord2[2, i]
        dz = coord1[3, i] - coord2[3, i]
        dist2 = dx^2.0 + dy^2.0 + dz^2.0
        total += dist2 * weight[i]
    end
    Σw = sum(weight)
    rms = sqrt(total / Σw)

    return rms
end

function default_pair(n::Int)

    ipair = zeros(Int, 2, n)
    for i = 1:n
        ipair[1, i] = i
        ipair[2, i] = i
    end

    return ipair
end

function default_weight(n::Int)

    weight = ones(dp, n)

    return weight
end

function coord_center!(coord1::Matrix{dp}, coord2::Matrix{dp}, ipair::Matrix{Int}, weight::Vector{dp})

    n = length(weight)
    xcent1 = 0.0
    ycent1 = 0.0
    zcent1 = 0.0
    xcent2 = 0.0
    ycent2 = 0.0
    zcent2 = 0.0
    for i = 1:n
        at1 = ipair[1, i]
        at2 = ipair[2, i]
        #coord1
        xcent1 += coord1[1, at1] * weight[i]
        ycent1 += coord1[2, at1] * weight[i]
        zcent1 += coord1[3, at1] * weight[i]
        #coord2
        xcent2 += coord2[1, at2] * weight[i]
        ycent2 += coord2[2, at2] * weight[i]
        zcent2 += coord2[3, at2] * weight[i]
    end
    Σw = sum(weight)
    xcent1 /= Σw
    ycent1 /= Σw
    zcent1 /= Σw
    xcent2 /= Σw
    ycent2 /= Σw
    zcent2 /= Σw
    for i = 1:n
        coord1[1, i] -= xcent1
        coord1[2, i] -= ycent1
        coord1[3, i] -= zcent1
        coord2[1, i] -= xcent2
        coord2[2, i] -= ycent2
        coord2[3, i] -= zcent2
    end

    return coord1, coord2
end

function symm_matrix(coord1::Matrix{dp}, coord2::Matrix{dp}, ipair::Matrix{Int}, weight::Vector{dp})

    n = length(weight)
    xp = zeros(dp, n)
    xm = zeros(dp, n)
    yp = zeros(dp, n)
    ym = zeros(dp, n)
    zp = zeros(dp, n)
    zm = zeros(dp, n)
    for i = 1:n
        at1 = ipair[1, i]
        at2 = ipair[2, i]
        xp[i] = (coord2[1, at2] + coord1[1, at1]) * weight[i]
        xm[i] = (coord2[1, at2] - coord1[1, at1]) * weight[i]
        yp[i] = (coord2[2, at2] + coord1[2, at1]) * weight[i]
        ym[i] = (coord2[2, at2] - coord1[2, at1]) * weight[i]
        zp[i] = (coord2[3, at2] + coord1[3, at1]) * weight[i]
        zm[i] = (coord2[3, at2] - coord1[3, at1]) * weight[i]
    end
    m = zeros(dp, 4, 4)
    for i = 1:n
        m[1, 1] += xm[i]^2.0 + ym[i]^2.0 + zm[i]^2.0
        m[1, 2] += yp[i] * zm[i] - ym[i] * zp[i]
        m[1, 3] += xm[i] * zp[i] - xp[i] * zm[i]
        m[1, 4] += xp[i] * ym[i] - xm[i] * yp[i]
        m[2, 2] += yp[i]^2.0 + zp[i]^2.0 + xm[i]^2.0
        m[2, 3] += xm[i] * ym[i] - xp[i] * yp[i]
        m[2, 4] += xm[i] * zm[i] - xp[i] * zp[i]
        m[3, 3] += xp[i]^2.0 + zp[i]^2.0 + ym[i]^2.0
        m[3, 4] += ym[i] * zm[i] - yp[i] * zp[i]
        m[4, 4] += xp[i]^2.0 + yp[i]^2.0 + zm[i]^2.0
    end
    # we need the full symmetric matrix
    m[2, 1] = m[1, 2]
    m[3, 2] = m[2, 3]
    m[3, 1] = m[1, 3]
    m[4, 3] = m[3, 4]
    m[4, 2] = m[2, 4]
    m[4, 1] = m[1, 4]

    return m
end

function rot_matrix(q::Vector{dp})

    rot = zeros(dp, 3, 3)

    rot[1,1] = q[1]^2 + q[2]^2 - q[3]^2 - q[4]^2
    rot[2,1] = 2.0 * (q[2] * q[3] - q[1] * q[4])
    rot[3,1] = 2.0 * (q[2] * q[4] + q[1] * q[3])

    rot[1,2] = 2.0 * (q[2] * q[3] + q[1] * q[4])
    rot[2,2] = q[1]^2 + q[3]^2 - q[2]^2 - q[4]^2
    rot[3,2] = 2.0 * (q[3] * q[4] - q[1] * q[2])

    rot[1,3] = 2.0 * (q[2] * q[4] - q[1] * q[3])
    rot[2,3] = 2.0 * (q[3] * q[4] + q[1] * q[2])
    rot[3,3] = q[1]^2 + q[4]^2 - q[2]^2 - q[3]^2

    return rot
end
