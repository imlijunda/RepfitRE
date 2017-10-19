function read_xyz(filepath::AbstractString)

    f = open(filepath, "r")
    lines = readlines(f)
    n = length(lines)
    close(f)

    #1 n_atom
    tmp = split(lines[1])
    n_atom = parse(Int, tmp[1])
    pbc = false

    #2 comment
    comment = strip(lines[2])

    #3-n coordinates
    coord = zeros(dp, 3, n_atom)
    sp_list = zeros(Int, n_atom)
    idx = 3
    for i = 1:n_atom
        s = split(lines[idx])
        sp_list[i] = atomic_num(s[1])
        coord[:, i] = parse_list(dp, @view(s[2:4]))
        idx += 1
    end

    return Geometry(coord, sp_list, comment)
end

function write_xyz(filepath::AbstractString, geo::Geometry)

    coord = geo.lattice * geo.coord
    b = 0

    open(filepath, "w") do f

        #1 n_atom
        n_atom = size(coord, 2)
        s = @sprintf("%d\n", n_atom)
        b += write(f, s)

        #2 comment
        s = @sprintf("%s\n", geo.comment)
        b += write(f, s)

        #3-n coordinates
        for i = 1:n_atom
            sym = atomic_sym(geo.sp[i])
            s = @sprintf("%s  %.12f  %.12f  %.12f\n", sym, coord[:, i]...)
            b += write(f, s)
        end
    end

    return b
end

function read_poscar(filepath::AbstractString)

    f = open(filepath, "r")
    lines = readlines(f)
    n = length(lines)
    close(f)

    #1 comment
    comment = strip(lines[1])

    #2 scaling
    scaling = parse(dp, lines[2])

    #3-5 lattice vector
    lattice = zeros(dp, 3, 3)
    for i = 3:5
        tmp = split(lines[i])
        lattice[:, i-2] = parse_list(dp, tmp) .* scaling
    end

    #6/7 number of atoms per atomic species or species names
    idx = 6
    tmp = split(lines[idx])
    n_species = length(tmp)
    sp_symbol = String[]
    sp_number = zeros(Int, n_species)
    sp_count = zeros(Int, n_species)
    if isnull(tryparse(Int, tmp[1]))
        #species names are specified
        sp_symbol = split(lines[idx])
        sp_number = atomic_num(sp_symbol)
        idx += 1
        tmp = split(lines[idx])
    end
    sp_count = parse_list(Int, tmp)
    n_atom = sum(sp_count)

    #fill sp, if not specified in POSCAR, zeros will be filled
    sp = zeros(Int, n_atom)
    ptr = 1
    for i = 1:n_species
        rng = ptr : ptr + sp_count[i] - 1
        sp[rng] = sp_number[i]
        ptr += sp_count[i]
    end

    #7/8 Selective or Cartesian or Direct
    idx += 1
    tmp = split(lines[idx])[1]
    # Constraints are not parse since technically they are not geometric information
    #selective = false
    if (tmp[1] == 'S') || (tmp[1] == 's')
        #selective = true
        idx += 1
    end
    #Cartesian or Direct
    tmp = split(lines[idx])[1]
    cartesian = (tmp[1] == 'C') || (tmp[1] == 'c')
    idx += 1

    #Coordinates
    coord = zeros(dp, 3, n_atom)
    for i = 1:n_atom
        tmp = split(lines[idx])
        coord[:, i] = parse_list(dp, @view(tmp[1:3]))
        idx += 1
    end
    #Velocities are ignored
    if cartesian
        #Convert to fractional coordinates
        coord = lattice \ coord
    end

    return Geometry(lattice, coord, sp, comment)
end

function write_poscar(filepath::AbstractString, geo::Geometry, constraints::Matrix{Bool} = Matrix{Bool}(0,0))

    b = 0
    open(filepath, "w") do f
        #comment and scaling (always 1.0)
        line = @sprintf("%s \n 1.0 \n", geo.comment)
        b += write(f, line)

        #lattice vector
        for i = 1:3
            line = @sprintf("%.12f  %.12f  %.12f \n", geo.lattice[:, i]...)
            b += write(f, line)
        end

        #now figure out sp_symbol and sp_count
        sp_symbol = String[]
        sp_count = Int[]
        push!(sp_symbol, atomic_sym(geo.sp[1]))
        push!(sp_count, 1)
        n_species = 1
        n_atom = size(geo.coord, 2)
        for i = 2:n_atom
            if geo.sp[i] != geo.sp[i-1]
                #new species
                push!(sp_symbol, atomic_sym(geo.sp[i]))
                push!(sp_count, 1)
                n_species += 1
            else
                #same species as previous
                sp_count[end] += 1
            end
        end

        #write species symbols
        for i = 1:n_species
            line = @sprintf("  %s", sp_symbol[i])
            b += write(f, line)
        end
        b += write(f, "\n")
        #write species count
        for i = 1:n_species
            line = @sprintf("  %d", sp_count[i])
            b += write(f, line)
        end
        b += write(f, "\n")

        selective = !isempty(constraints)
        if selective
            b += write(f, "Selective dynamics\n")
        end

        #write coordinates
        b += write(f, "Direct\n")
        if !selective
            for i = 1:n_atom
                line = @sprintf("%.12f  %.12f  %.12f\n", geo.coord[:, i]...)
                b += write(f, line)
            end
        else
            for i = 1:n_atom
                line = @sprintf("%.12f  %.12f  %.12f ", geo.coord[:, i]...)
                b += write(f, line)
                for j = 1:3
                    b += write(f, constraints[j, i] ? " F":" T")
                end
                b += write(f, "\n")
            end
        end
    end

    return b
end

function read_gen(filepath::AbstractString)

    f = open(filepath, "r")
    lines = readlines(f)
    n = length(lines)
    close(f)

    comment = ""

    #1 n_atom and F/C
    tmp = split(lines[1])
    n_atom = parse(Int, tmp[1])
    pbc = (tmp[2] == "F") || (tmp[2] == "f")

    #2 species list
    sp_symbol = split(lines[2])
    sp_number = atomic_num(sp_symbol)

    coord = zeros(dp, 3, n_atom)
    sp = zeros(Int, n_atom)
    idx = 3
    for i = 1:n_atom
        tmp = split(lines[idx])
        sp[i] = sp_number[parse(Int, tmp[2])]
        coord[:, i] = parse_list(dp, @view(tmp[3:5]))
        idx += 1
    end

    lattice = zeros(dp, 3, 3)
    if pbc
        #ignore origin shift
        idx += 1
        #lattic vector
        for i = 1:3
            tmp = split(lines[idx])
            lattice[:, i] = parse_list(dp, tmp)
            idx += 1
        end
        reslt = Geometry(lattice, coord, sp, comment)
    else
        reslt = Geometry(coord, sp, comment)
    end

    return reslt
end

function write_gen(filepath::AbstractString, geo::Geometry)

    b = 0
    open(filepath, "w") do f

        # number of atoms and structure type
        n_atom = size(geo.coord, 2)
        if geo.pbc
            line = @sprintf("%d F\n", n_atom)
        else
            line = @sprintf("%d C\n", n_atom)
        end
        b += write(f, line)

        #now figure out sp_symbol and sp_number
        sp_symbol = String[]
        sp_number = Int[]
        n_species = 1
        push!(sp_number, geo.sp[1])
        for i = 2:n_atom
            if locate(geo.sp[i], sp_number) == 0
                push!(sp_number, geo.sp[i])
                n_species += 1
            end
        end
        sp_symbol = atomic_sym(sp_number);
        #write species symbols
        for i = 1:n_species
            line = @sprintf("  %s", sp_symbol[i])
            b += write(f, line)
        end
        b += write(f, "\n")

        #write coordinates
        for i = 1:n_atom
            sp_idx = locate(geo.sp[i], sp_number)
            line = @sprintf("%d  %d  %.12f  %.12f  %.12f\n", i, sp_idx, geo.coord[:, i]...)
            b += write(f, line)
        end

        if geo.pbc
            #write lattice vectors
            b += write(f, "0.0  0.0  0.0\n")
            for i = 1:3
                line = @sprintf("%.12f  %.12f  %.12f\n", geo.lattice[:, i]...)
                b += write(f, line)
            end
        end

        #write comment
        line = @sprintf("\n#%s\n", geo.comment)
        b += write(f, line)
    end

    return b
end
