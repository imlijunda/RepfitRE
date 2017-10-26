const SCALE_PREFIX = "scaled-"

function write_structures(wd::AbstractString, geo::Geometry; vasp = "POSCAR", dftb = "geo.gen")

    b = 0
    fname = wd_fname(wd, vasp)
    b += write_poscar(fname, geo)
    fname = wd_fname(wd, dftb)
    b += write_gen(fname, geo)

    return b
end

function write_repscan_structures(wd::AbstractString, rsc::RepScale; vasp = "POSCAR", dftb = "geo.gen")

    # write root structure first
    mkpath(wd)
    write_structures(wd, rsc.rep.nei.geo; vasp = vasp, dftb = dftb)

    # files written
    f = 2
    # write scaled structures
    for i = 1:length(rsc.scaled)
        new_wd = wd_tagfname(wd, i, SCALE_PREFIX)
        new_rsc = rsc.scaled[i]
        f += write_repscan_structures(new_wd, new_rsc; vasp = vasp, dftb = dftb)
    end

    return f
end

function read_repscan_energies!(wd::AbstractString, rsc::RepScale; vasp = "OUTCAR", dftb = "detailed.out")

    # read root energies
    fname = wd_fname(wd, vasp)
    rsc.edft = vasp_esigma0(fname)

    fname = wd_fname(wd, dftb)
    dftb_tot = dftb_mermin(fname)
    dftb_rep = dftb_repulsive(fname)
    rsc.eband = dftb_tot - dftb_rep

    rsc.tot_erep = rsc.edft - rsc.eband

    # read scaled structures
    for i = 1:length(rsc.scaled)
        new_wd = wd_tagfname(wd, i, SCALE_PREFIX)
        read_repscan_energies!(new_wd, rsc.scaled[i]; vasp = vasp, dftb = dftb)
    end

    return rsc
end
