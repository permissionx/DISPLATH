# =============================================================================
# Example: loading an arbitrary structure from a file (a C60 fullerene)
#          and irradiating it with a single Ne ion  (static mode)
# -----------------------------------------------------------------------------
# Instead of generating a crystal from lattice vectors, the target is read from
# a LAMMPS data file ("C60.data", atomic style). This is the route to use for
# any custom / non-periodic / experimentally derived geometry (molecules,
# clusters, defected supercells, supported stacks, ...).
#
# Loading from a file requires STATIC mode (IS_DYNAMIC_LOAD = false).
#
# Run with:   julia main.jl
# =============================================================================

const IS_DYNAMIC_LOAD = false             # required for structure-from-file
home = ENV["ARCS_HOME"]
include(home * "/src/DISPLATH.jl")

seed = 43
const THREAD_RNG = [StableRNG(seed + t) for t in 1:Threads.nthreads()]

# --- parameters ---------------------------------------------------------------
pMax = 3.5
vacancyRecoverDistance = 0.0
parameters = Parameters(pMax, vacancyRecoverDistance;
                        isDumpInCascade = true,   # write a frame per collision
                        isNonQnl        = true)

# --- target read from file ----------------------------------------------------
typeDict = Dict(
    1 => Element("C",  22.0, 7.9),   # (Ed, binding energy) in eV
    2 => Element("Ne", 0.1,  0.1),   # projectile
)
inputGridVectors = [4.0 0.0 0.0; 0.0 4.0 0.0; 0.0 0.0 4.0]

# `replicate` tiles the cell along x/y/z (use [1,1,1] for none)
material  = Material("C60.data", typeDict, inputGridVectors, parameters;
                     replicate = [1, 1, 1])
simulator = Simulator(material, parameters)
Save!(simulator)
@dump "C60.dump" simulator.atoms          # initial structure

# --- single 1 keV Ne impact, aimed at the cage from above ---------------------
energy      = 1000.0
ionPosition = RandomPointInCircle(5.0) + [10.2306, 10.4537, 15.8964]
Irradiation!(simulator, energy, ionPosition, [0.0, 0.0, -1.0], 2, parameters)
@dump "C60_after.dump" simulator.atoms     # structure after the cascade
