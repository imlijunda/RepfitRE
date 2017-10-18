# 2D surface scan

# 2D scanning surface
# A surface is simply defined by vector a and vector b
mutable struct Surface

    origin::Vector{dp}
    vecta::Vector{dp}
    vectb::Vector{dp}
    sh::Vector{dp}

    function Surface(geo::Geometry, at0::Int, at1::Int, at2::Int)

        s = new()
        s.origin = geo.coord[:, at0]
        s.vecta = geo.coord[:, at1] - s.origin
        s.vectb = geo.coord[:, at2] - s.origin
        s.sh = zeros(dp, 3)

        return s
    end
end

# Verticle shift of the surface, positive direction is defined by right hand
# rule of vect a-b
function shift_surface!(s::Surface, distance::dp)

    n = cross(s.vecta, s.vectb)
    n = n / norm(n)
    s.sh = n * distance

    return s
end

mutable struct SurfaceScan

    surf::Surface

    na::Int
    nb::Int
    delta::dp

    spfunc::Function

    function SurfaceScan(surf::Surface, na::Int, nb::Int; delta::dp = 0.0, spfunc = linspace)

        ss = new()
        ss.surf = surf
        ss.na = na
        ss.nb = nb
        ss.delta = 0.0
        ss.spfunc = spfunc

        return ss
    end
end

function coord_mesh(ss::SurfaceScan)

    va = ss.spfunc(0.0, 1.0, ss.na)
    va += ss.delta
    vb = ss.spfunc(0.0, 1.0, ss.nb)
    vb += ss.delta

    cmesh = Matrix{Vector{dp}}(ss.na, ss.nb)
    (xx, yy) = ndgrid(va, vb)
    for i = 1:ss.na
        for j = 1:ss.nb
            cmesh[i, j] = ss.surf.vecta * xx[i, j] + ss.surf.vectb * yy[i, j] + ss.surf.sh
        end
    end

    return cmesh
end

function geo_mesh(ss::SurfaceScan, geo::Geometry, adsorbate::AtomList)

    cmesh = coord_mesh(ss)
    gmesh = Matrix{Geometry}(ss.na, ss.nb)

    for i = 1:ss.na
        for j = 1:ss.nb
            newgeo = Geometry(geo)
            move!(newgeo, adsorbate, cmesh[i, j])
            gmesh[i, j] = newgeo
        end
    end

    return gmesh
end
