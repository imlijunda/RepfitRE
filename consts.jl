# Type alias
const dp = Float64

# ~machine epsilon of 0.1 as single precision fp
const PRECIS = 1e-8

# values taken from NIST Reference website: 2014 CODATA recommended values
# http://physics.nist.gov/cuu/Constants/index.html
const ang_per_bohr = 0.52917721067
const bohr_per_ang = 1.0 / ang_per_bohr
const ev_per_hartree = 27.21138602
const hartree_per_ev = 1.0 / ev_per_hartree
# Boltzmann const
const ev_per_kelvin = 8.617330350e-5
const kelvin_per_ev = 1.0 / ev_per_kelvin
# Periodic Table
const p_table = ["H","He","Li","Be","B","C","N","O","F","Ne","Na","Mg","Al","Si",
        "P","S","Cl","K","Ar","Ca","Sc","Ti","V","Cr","Mn","Fe","Ni","Co","Cu",
        "Zn","Ga","Ge","As","Se","Br","Kr","Rb","Sr","Y","Zr","Nb","Mo","Tc",
        "Ru","Rh","Pd","Ag","Cd","In","Sn","Sb","I","Te","Xe","Cs","Ba","La",
        "Ce","Pr","Nd","Pm","Sm","Eu","Gd","Tb","Dy","Ho","Er","Tm","Yb","Lu",
        "Hf","Ta","W","Re","Os","Ir","Pt","Au","Hg","Tl","Pb","Bi","Po","At",
        "Rn","Fr","Ra","Ac","Pa","Th","Np","U","Am","Pu","Cm","Bk","Cf","Es",
        "Fm","Md","No","Lr","Rf","dp","Sg","Hs","Bh","Mt","Ds","Rg","Cn","Uut",
        "Fl","Uup","Lv","Uus","Uuo"]
