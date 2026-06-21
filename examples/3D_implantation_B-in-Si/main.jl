# =============================================================================
# Example: B implantation into 3D crystalline Si(100)  (dynamic-loading mode)
# -----------------------------------------------------------------------------
# Reproduces the type of calculation shown in Fig. 3-4 of the ASIM paper: the
# depth distribution R_p of boron implanted into Si is collected over many
# incident ions and written to "R_p.csv".
#
# Because a keV/MeV implant must cover a micrometer-scale penetration range
# (>10^7 atoms), the global lattice is NOT instantiated up front. Instead the
# DYNAMIC-LOADING mode is used: atoms are created on the fly only inside the
# cells an active atom actually visits, and are periodically released, while
# the produced defects are retained.  The mode is selected by the global flag
# `IS_DYNAMIC_LOAD`, which MUST be defined before `DISPLATH.jl` is included.
#
# Run with:   julia main.jl              (ARCS_HOME / ARCS_REPO must be set)
# =============================================================================

const IS_DYNAMIC_LOAD = true              # <-- enable dynamic loading
home = ENV["ARCS_HOME"]
include(home * "/src/DISPLATH.jl")

# Reproducible per-thread random streams (required by ASIM)
seed = 43
const THREAD_RNG = [StableRNG(seed + t) for t in 1:Threads.nthreads()]

# --- number of incident ions and target lateral size from the target fluence -
NI   = 10000                              # number of incident ions
flux = 5e14 * 1e-16                       # fluence (ions / Angstrom^2)

# --- collision / physics parameters ------------------------------------------
a = 5.431                                 # Si lattice constant (Angstrom)
pMax = a / 3                              # impact-parameter cut-off
vacancyRecoverDistance = 0.0             # no spontaneous recombination here

# Note: typeDict is supplied to `Material` (below), not to `Parameters`.
typeDict = Dict(
    1 => Element("Si", 20.0, 10.0),       # target: (Ed, binding energy) in eV
    2 => Element("B",   0.1,  0.1),       # projectile
)

parameters = Parameters(pMax, vacancyRecoverDistance;
                        temperature       = 300.0,    # K, thermal vibrations on
                        DebyeTemperature  = 519.0,    # K, Debye model for Si
                        stopEnergy        = 10.0,     # eV, transport cut-off
                        nCascadeEveryLoad = 1000,     # cascades between releases
                        isAmorphous       = false,    # crystalline target
                        amorphousLength   = 5.0,      # thin amorphous top layer
                        maxRSS            = 100)      # memory budget (GB)

# --- target lattice: Si conventional diamond cell (8 atoms) ------------------
boxL = ceil(Int, sqrt(NI / flux / a / a))   # lateral cells to match the fluence
boxL = max(boxL, 3)                          # the link-cell grid needs >= 3 cells
log_info("boxL: $boxL")
primaryVectors = [a 0.0 0.0; 0.0 a 0.0; 0.0 0.0 a]
boxSizes       = [boxL, boxL, 15005]
latticeRanges  = [0 boxL; 0 boxL; 2 15000]
basis = [0.0  0.0  0.0;
         0.5  0.5  0.0;
         0.5  0.0  0.5;
         0.0  0.5  0.5;
         0.25 0.25 0.25;
         0.75 0.75 0.25;
         0.75 0.25 0.75;
         0.25 0.75 0.75]
basisTypes = [1, 1, 1, 1, 1, 1, 1, 1]       # all Si
# In dynamic mode the grid spacing must be an integer multiple of the unit cell.
inputGridVectors = [a 0.0 0.0; 0.0 a 0.0; 0.0 0.0 a]

material  = Material(primaryVectors, latticeRanges, basisTypes, basis, typeDict,
                     boxSizes, inputGridVectors, parameters)
simulator = Simulator(material, parameters)
Restore!(simulator)

# --- irradiation conditions ---------------------------------------------------
energy = 10_000.0                          # 10 keV
phi    = 30 / 180 * pi                      # rotation angle
theta  = 7  / 180 * pi                      # tilt angle (0deg = channelling)
v0     = [cos(phi) * sin(theta), sin(phi) * sin(theta), -cos(theta)]
zTop   = latticeRanges[3, 2] * a - 2       # start just below the surface

for i in 1:NI
    @show i
    rp          = RandomInSquare(boxSizes[1] * a, boxSizes[2] * a)   # random (x,y)
    ionPosition = Vector{Float64}(rp) + [0.0, 0.0, zTop]
    v           = RandomlyDeviatedVector(v0, 0.5 / 180 * pi)          # 0.5deg beam divergence
    ion = Atom(2, ionPosition, parameters)
    SetVelocityDirection!(ion, v)
    SetEnergy!(ion, energy)
    push!(simulator, ion)
    @time Cascade!(ion, simulator)

    # projected range = distance travelled below the surface
    @record "R_p.csv" "$(zTop - ion.coordinate[3] + 2)"
end
