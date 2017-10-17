const AtomList = Vector{Int}

mutable struct Geometry

    # periodic boundary condition
    pbc::Bool

    # lattice vectors and coordinates are all stored by column vectors
    # lattice vector, should keep as I marix if not matrix for consistency
    lattice::Matrix{dp}
    # coordinates
    coord::Matrix{dp}
    # atomic species corresponding to coord
    sp::Vector{Int}

    comment::String
    tag::Int

    # use constructor to copy
    function Geometry(geo::Geometry)

        newgeo = new()

        newgeo.pbc = geo.pbc
        newgeo.lattice = copy(geo.lattice)
        newgeo.coord = copy(geo.coord)
        newgeo.sp = copy(geo.sp)
        newgeo.comment = geo.comment
        newgeo.tag = geo.tag

        return newgeo
    end

    # create new geometry from selected atom list
    function Geometry(geo::Geometry, al::AtomList, comment::AbstractString, tag::Int = 0)

        newgeo = new()

        newgeo.pbc = geo.pbc
        newgeo.lattice = copy(geo.lattice)
        n = length(al)
        newgeo.coord = zeros(dp, 3, n)
        newgeo.sp = zeros(dp, n)
        for i = 1:n
            for j = 1:3
                newgeo.coord[j, i] = geo.coord[j, al[i]]
            end
            newgeo.sp[i] = geo.sp[al[i]]
        end
        newgeo.comment = string(comment)
        newgeo.tag = tag

        return newgeo
    end

    # create new periodic boundary condition geometry
    function Geometry(lattice::Matrix{dp}, coord::Matrix{dp}, sp::Vector{Int}, comment::AbstractString, tag::Int = 0)

        geo = new()

        geo.pbc = true
        geo.lattice = zeros(dp, 3, 3)
        geo.coord = zeros(dp, size(coord))
        n = size(coord, 2)
        for i = 1:3
            for j = 1:3
                geo.lattice[i, j] = lattice[i, j]
            end
            for j = 1:n
                geo.coord[i, j] = coord[i, j]
            end
        end
        geo.sp = copy(sp)
        geo.comment = string(comment)
        geo.tag = tag

        return geo
    end

    # create new non-periodic geometry
    function Geometry(coord::Matrix{dp}, sp::Vector{Int}, comment::AbstractString, tag::Int = 0)

        geo = new()

        geo.pbc = false
        geo.lattice = eye(dp, 3)
        geo.coord = zeros(dp, size(coord))
        n = size(coord, 2)
        for i = 1:n
            for j = 1:3
                geo.coord[j, i] = coord[j, i]
            end
        end
        geo.sp = copy(sp)
        geo.comment = string(comment)
        geo.tag = tag

        return geo
    end
end

# Geometry functions

include("quaternion.jl")
function goemtry_rms(geo1::Geometry, geo2::Geometry)

    if geo1.pbc
        coord1 = goe1.lattice * goe1.coord
    else
        coord1 = copy(geo1.coord)
    end
    if geo2.pbc
        coord2 = geo2.lattice * geo2.coord
    else
        coord2 = copy(geo2.coord)
    end
    n1 = size(coord1, 2)
    n2 = size(coord2, 2)

    if n1 != n2
        warn("Comparing two geometries with different numbers of atoms, ignoring extra atoms")
        n = min(n1, n2)
    else
        n = n1
    end
    ipair = default_pair(n)
    weight = default_weight(n)

    quaternion!(coord1, coord2, ipair, weight)
    return coord_rms(coord1, coord2, ipair, weight)
end

function select_atoms_species(geo::Geometry, sp::Union{Int, Vector{Int}})

    atom_list = Int[]
    n = size(geo.coord, 2)
    for i = 1:n
        if in(geo.sp[i], sp)
            push!(atom_list, i)
        end
    end

    return atom_list
end

# Calculate the component of projected vector on the main vector direction
scalar_projection(main::Vector{dp}, projected::Vector{dp}) = dot(main, projected) / norm(main)
function select_atoms_axis(geo::Geometry, v::Int, threshold::dp, sp::Union{Int, Vector{Int}} = Int[])

    ax = geo.lattice[:, v]
    len_ax = norm(ax)
    emptysp = isempty(sp)
    n = size(geo.coord, 2)

    atom_list = Int[]
    tmp = zeros(dp, 3)

    for i = 1:n
        if emptysp || in(geo.sp[i], sp)
            #copy coordinates to tmp so we don't create new local variable in every loop
            for j = 1:3
                tmp[j] = geo.coord[j, i]
            end
            comp_percentage = scalar_projection(ax, tmp) / len_ax
            if comp_percentage > threshold
                push!(atom_list, i)
            end
        end
    end

    return atom_list
end
function select_atoms_axis(geo::Geometry, v::Vector{dp}, threshold::dp, sp::Union{Int, Vector{Int}} = Int[])

    len_v = norm(v)
    emptysp = isempty(sp)
    n = size(geo.coord, 2)

    atom_list = Int[]
    tmp = zeros(dp, 3)

    for i = 1:n
        if emptysp || in(geo.sp[i], sp)
            #copy coordinates to tmp so we don't create new local variable in every loop
            for j = 1:3
                tmp[j] = geo.coord[j, i]
            end
            comp_percentage = scalar_projection(v, tmp) / len_v
            if comp_percentage >= threshold
                push!(atom_list, i)
            end
        end
    end

    return atom_list
end

function select_atoms_sphere(geo::Geometry, p::Vector{dp}, r::dp, sp::Union{Int, Vector{Int}} = Int[])

    #get absolute coordinates
    coord = geo.lattice * geo.coord
    emptysp = isempty(sp)
    n = size(coord, 2)

    atom_list = Int[]
    tmp = zeros(dp, 3)

    for i = 1:n
        if emptysp || in(geo.sp[i], sp)
            #copy coordinates to tmp so we don't create new local variable in every loop
            for j = 1:3
                #not using dist() to avoid creating new local variable
                tmp[j] = coord[j, i] - p[j]
            end
            if norm(tmp) <= r
                push!(atom_list, i)
            end
        end
    end

    return atom_list
end
