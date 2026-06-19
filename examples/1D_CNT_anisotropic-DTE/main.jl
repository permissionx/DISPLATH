# =============================================================================
# Example: Ar irradiation of a 1D (10,10) single-walled carbon nanotube
#          with a USER-DEFINED, direction-dependent displacement threshold Ed
# -----------------------------------------------------------------------------
# This reproduces the CNT case of the ASIM paper (Fig. 8). The displacement
# threshold of a recoiling C atom depends on whether it is knocked OUTWARD
# (away from the tube axis, Ed = 14 eV) or INWARD (towards the axis, Ed = 25 eV)
# [Merrill et al., PRB 92, 075404 (2015)].
#
# Such a rule cannot be expressed by a single per-species number. ASIM exposes
# it through `DTEMode = 3` ("custom"): the solver then calls two user-supplied
# functions, `GetDTECustom` and `GetBDECustom`, for every potential recoil.
# Define them at top level in this script -- no change to the core code needed.
#
# The structure is read from a LAMMPS data file; STATIC loading is used.
#
# Run with:   julia main.jl
# =============================================================================

const IS_DYNAMIC_LOAD = false             # structure-from-file => static mode
home = ENV["ARCS_HOME"]
include(home * "/src/DISPLATH.jl")

seed = 43
const THREAD_RNG = [StableRNG(seed + t) for t in 1:Threads.nthreads()]

# -----------------------------------------------------------------------------
# Custom displacement-threshold rule (called by the solver because DTEMode == 3)
# -----------------------------------------------------------------------------
# The tube axis is along x; the cross-section lies in the (y, z) plane. A recoil
# moves OUTWARD when the radial vector r (from the axis to the atom) and the
# in-plane velocity v point the same way, i.e. dot(r, v) > 0.
# `TUBE_CENTER` is filled in after the structure is loaded (see below).
const TUBE_CENTER = Float64[0.0, 0.0]

function GetDTECustom(atom::Atom, simulator::Simulator)
    r = atom.coordinate[2:3] .- TUBE_CENTER       # radial direction in (y,z)
    v = atom.velocityDirection[2:3]               # in-plane recoil direction
    return dot(r, v) > 0 ? 14.0 : 25.0            # outward : inward  (eV)
end

GetBDECustom(atom::Atom, simulator::Simulator) = atom.bde

# -----------------------------------------------------------------------------
# Parameters and target
# -----------------------------------------------------------------------------
pMax = 1.2
vacancyRecoverDistance = 2.88
parameters = Parameters(pMax, vacancyRecoverDistance;
                        DTEMode  = 3,        # <-- use GetDTECustom / GetBDECustom
                        isNonQnl = true)     # local electronic stopping only

typeDict = Dict(
    1 => Element("C",  20.0, 7.9),    # (Ed placeholder; overridden by DTEMode 3)
    2 => Element("Ar", 20.0, 7.9),    # projectile
)
inputGridVectors = [4.0 0.0 0.0; 0.0 4.0 0.0; 0.0 0.0 4.0]

material  = Material("cnt_10_10.data", typeDict, inputGridVectors, parameters)
simulator = Simulator(material, parameters)
Save!(simulator)

# tube cross-section centre = mid-point of the box in (y, z)
TUBE_CENTER[1] = simulator.box.vectors[2, 2] / 2
TUBE_CENTER[2] = simulator.box.vectors[3, 3] / 2

@dump "cnt.dump" simulator.atoms

# -----------------------------------------------------------------------------
# Sputtering yield vs. Ar energy
# -----------------------------------------------------------------------------
function CountVacancies2D(simulator::Simulator)
    nV = 0
    for atomIndex in simulator.displacedAtoms
        simulator.atoms[atomIndex].isAlive || (nV += 1)
    end
    return nV
end

N = 2000                                   # incident ions per energy
xlen = simulator.box.vectors[1, 1]         # tube axis (periodic)
ylen = simulator.box.vectors[2, 2]
zlen = simulator.box.vectors[3, 3]
for energy in 100.0:100.0:2000.0
    nV = 0.0
    for _ in 1:N
        Restore!(simulator)
        # uniform impact point over the tube footprint, ion travelling in -z
        rp = RandomInSquare(xlen, ylen)
        ionPosition = [rp[1], rp[2], zlen + 5.0]
        Irradiation!(simulator, energy, ionPosition, [0.0, 0.0, -1.0], 2, parameters)
        nV += CountVacancies2D(simulator)
    end
    nV /= N
    @show energy, nV
    @record "nV.csv" "$(energy),$(nV)" "energy,nV"
end
