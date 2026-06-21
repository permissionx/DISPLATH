# ASIM — Atomic-scale Simulation on Ion Irradiation of Matter

ASIM is an open-source, Julia-based **Binary Collision Approximation (BCA)**
framework for simulating ion irradiation of materials of *any* dimensionality —
from micrometre-scale 3D bulk (e.g. ion implantation into Si) to 2D layers
(graphene, h-BN), 1D nanostructures (carbon nanotubes) and 0D clusters
(fullerenes), as well as their combinations (e.g. substrate-supported 2D
sheets). It models both nuclear and electronic stopping, explicitly treats
simultaneous (multi-target) collisions, and outputs defect statistics, defect
morphologies and the full structural evolution of the target.

The distinguishing feature of ASIM is **programmability**: the irradiation
scenario, the output observables and — most importantly — the displacement
threshold rule can all be supplied as short Julia functions in the input file,
*without modifying the core code*. This README documents how to install ASIM,
how a simulation is built, and how to customise model systems, irradiation
conditions and displacement-threshold rules. Worked, runnable examples are
provided under [`examples/`](examples/).

---

## 1. Requirements

* **Julia ≥ 1.9** (developed and tested with the build under
  Julia 1.11.7).
* External Julia packages (the only ones that must be installed; everything
  else is from the standard library):

  | Package          | Used for                                  |
  | ---------------- | ----------------------------------------- |
  | `StaticArrays`   | fast fixed-size vectors in the hot loop   |
  | `StableRNGs`     | reproducible per-thread random numbers    |
  | `QuadGK`         | the scattering integral (θ tabulation)    |
  | `Interpolations` | cubic interpolation of the θ/τ tables     |
  | `ProgressMeter`  | progress bars (`@showprogress`)           |
  | `Distributions`  | thermal displacement / beam divergence    |

Install them once with:

```julia
julia> import Pkg
julia> Pkg.add(["StaticArrays", "StableRNGs", "QuadGK",
                "Interpolations", "ProgressMeter", "Distributions"])
```

---

## 2. Installation and environment variables

ASIM is used by `include`-ing `src/DISPLATH.jl` from an input script. It relies
on two environment variables:

| Variable    | Meaning                                                                  |
| ----------- | ------------------------------------------------------------------------ |
| `ARCS_HOME` | path to the **repository root** (the folder that contains `src/`)         |
| `ARCS_REPO` | path to a directory holding `thetatau_repository/` and `dte_repository/`  |

`ARCS_REPO` is where ASIM **caches the scattering-angle tables** (`*.thetatau`).
The first time a projectile–target pair is encountered the table is computed
from the ZBL potential and written there; subsequent runs reuse it.

A helper script, [`install_arcs.sh`](install_arcs.sh), creates
`$HOME/.arcs/{thetatau_repository,dte_repository}` and appends the environment
variables (and the Julia path) to your `~/.bashrc`. You may instead set them by
hand, e.g.:

```bash
export ARCS_HOME=/path/to/DISPLATH          # the repo root
export ARCS_REPO=$HOME/.arcs                 # writable cache directory
mkdir -p $ARCS_REPO/thetatau_repository $ARCS_REPO/dte_repository
```

A set of pre-computed `*.thetatau` tables is also shipped in
[`thetatau_repository/`](thetatau_repository/); copy them into
`$ARCS_REPO/thetatau_repository/` to skip the first-run tabulation.

---

## 3. Running a simulation

Every input script is a plain Julia file. Run it with:

```bash
julia examples/load_structure_from_file_C60/main.jl      # single thread
julia -t 8 examples/3D_implantation_B-in-Si/main.jl      # 8 threads
```

A minimal script has the following five parts, in order:

```julia
const IS_DYNAMIC_LOAD = false                 # 1. choose the loading mode (see §5)
include(ENV["ARCS_HOME"] * "/src/DISPLATH.jl")

seed = 42                                      # 2. reproducible RNG (required)
const THREAD_RNG = [StableRNG(seed + t) for t in 1:Threads.nthreads()]

parameters = Parameters(pMax, vacancyRecoverDistance; ...)   # 3. global options
material   = Material(...)                                    # 4. build the target
simulator  = Simulator(material, parameters)

for i in 1:N                                   # 5. the irradiation scenario
    Restore!(simulator)
    Irradiation!(simulator, energy, position, direction, ionType, parameters)
    # ... collect observables ...
end
```

> **Note** `IS_DYNAMIC_LOAD` and `THREAD_RNG` are `const` globals that the core
> code reads; they **must** be defined as shown (the former *before* the
> `include`).

---

## 4. The API

### 4.1 Elements

```julia
Element(name::String, Ed::Float64, Eb::Float64)
```

`name` is the chemical symbol (e.g. `"Si"`, `"Ar"`); `Ed` is the displacement
threshold energy and `Eb` the binding energy, both in eV. Atomic mass, radius,
charge and the electronic-stopping parameters (α, β) are filled in from a
built-in table ([`src/elements.jl`](src/elements.jl)).

Elements are collected in a `typeDict` that maps an integer **type id** to an
`Element`. By convention the target species come first and the projectile last:

```julia
typeDict = Dict(
    1 => Element("Si", 20.0, 10.0),   # target
    2 => Element("B",  0.1,  0.1),    # projectile
)
```

### 4.2 Parameters — global options

```julia
parameters = Parameters(pMax, vacancyRecoverDistance; kwargs...)
```

* `pMax` — impact-parameter cut-off (Å); collisions with a larger impact
  parameter are ignored. A good default is ≈ half the lattice constant.
* `vacancyRecoverDistance` — an interstitial closer than this to a vacancy
  recombines (Å); `0.0` disables spontaneous recovery.

Commonly used keyword options:

| Keyword             | Default             | Meaning                                                   |
| ------------------- | ------------------- | --------------------------------------------------------- |
| `stopEnergy`        | `0.1`               | transport cut-off; an atom below this energy stops (eV)    |
| `DTEMode`           | `1`                 | displacement-threshold scheme — see §6                     |
| `DTEFile`           | `""`                | environment→Ed table for `DTEMode = 2`                     |
| `temperature`       | `0.0`               | target temperature for thermal vibrations (K)             |
| `DebyeTemperature`  | `519.0`             | Debye temperature for the vibration model (K)             |
| `isNonQnl`          | `false`             | use only local electronic stopping (no non-local term)    |
| `periodic`          | `[true,true,false]` | periodic boundaries along x/y/z                            |
| `isDumpInCascade`   | `false`             | write a dump frame after every collision (visualisation)  |
| `isAmorphous`       | `false`             | randomise lattice positions (amorphous target)            |
| `amorphousLength`   | —                   | thickness of an amorphous top layer (Å)                   |
| `nCascadeEveryLoad` | `100`               | (dynamic mode) cascades between memory releases           |
| `maxRSS`            | `20`                | (dynamic mode) resident-memory budget (GB)                |

> `typeDict` is **not** passed to `Parameters`; it is supplied to `Material`
> (below), which records it into the parameters object.

### 4.3 Material — building the target

Two constructors are available.

**(a) From lattice vectors** (ideal crystal):

```julia
material = Material(primaryVectors, latticeRanges, basisTypes, basis, typeDict,
                    boxSizes, inputGridVectors, parameters)
```

* `primaryVectors` — 3×3 unit-cell vectors (Å).
* `basis`, `basisTypes` — fractional positions of the basis atoms and their type ids.
* `latticeRanges` — integer `[lo hi]` per axis: which unit cells to populate.
* `boxSizes` — number of cells per axis defining the simulation box.
* `inputGridVectors` — 3×3 link-cell grid spacing (Å); must differ from
  `primaryVectors`.

**(b) From a structure file** (any custom geometry; **static mode only**):

```julia
material = Material(fileName, typeDict, inputGridVectors, parameters;
                    replicate = [1, 1, 1])
```

`fileName` is a LAMMPS data file (atomic style; see
`examples/load_structure_from_file_C60/C60.data`). `replicate` tiles the cell
along x/y/z.

Then:

```julia
simulator = Simulator(material, parameters)
```

### 4.4 Irradiating

The convenience driver fires one ion and runs its full cascade:

```julia
Irradiation!(simulator, energy, ionPosition, direction, ionType, parameters)
```

with `energy` in eV, `ionPosition`/`direction` as 3-vectors and `ionType` the
type id of the projectile. Equivalently, the low-level steps are:

```julia
ion = Atom(ionType, ionPosition, parameters)
SetVelocityDirection!(ion, direction)
SetEnergy!(ion, energy)
push!(simulator, ion)
Cascade!(ion, simulator)
```

Helpers for sampling impact points: `RandomInSquare(a, b)`,
`RandomPointInCircle(r)`, and `RandomlyDeviatedVector(v, divergence)` for beam
divergence.

### 4.5 State management, statistics and output

| Call                             | Purpose                                                 |
| -------------------------------- | ------------------------------------------------------- |
| `Save!(simulator)`               | store the pristine state (static mode)                  |
| `Restore!(simulator)`            | reset to the stored state before the next ion           |
| `DefectStatics(simulator)`       | returns `(interstitials, vacancies)`                    |
| `CountVacancies(simulator)`      | number of vacancies                                     |
| `@dump "file.dump" atoms`        | append a LAMMPS-style snapshot (e.g. `simulator.atoms`) |
| `@record "file.csv" value title` | append a line to a CSV (header `title` written once)    |

Use `Save!`/`Restore!` to accumulate statistics over many *independent* ions
(2D/1D/0D); omit them for a *sequential/cumulative* irradiation where damage
builds up (3D implantation at finite fluence).

---

## 5. Static vs. dynamic loading

The global flag `IS_DYNAMIC_LOAD` (set before the `include`) selects how the
target is held in memory:

* **Static** (`false`) — all atoms are instantiated up front and kept in the
  link-cell grid. Use for low-dimensional targets and small/medium crystals.
  Required when loading a structure from a file.
* **Dynamic** (`true`) — the global lattice is stored only implicitly; atoms are
  created on the fly only inside the cells an active atom visits, and are
  periodically released (`nCascadeEveryLoad`, `maxRSS`) while the produced
  defects are retained. This makes micrometre-scale 3D targets (>10⁷ atoms)
  tractable in modest memory.

The *same* input API is used in both modes; only the flag changes.

---

## 6. Customisation (custom systems, conditions and Ed rules)

This is the central design goal of ASIM. Three things are routinely customised.

### 6.1 Custom model systems

Any geometry can be supplied through a LAMMPS data file and
`Material(fileName, ...)` (static mode) — molecules, clusters, defected
supercells, nanotubes, or a 2D sheet stacked on a substrate. See
[`examples/load_structure_from_file_C60`](examples/load_structure_from_file_C60).

### 6.2 Custom irradiation conditions

Because the scenario is just Julia, arbitrarily complex irradiation protocols
are expressed as ordinary control flow: scan energies/angles/fluences, mix ion
species, ramp the energy, sample impact positions from any distribution, and
record any observable. See the energy scans in the examples.

### 6.3 Custom displacement-threshold rules (`DTEMode`)

ASIM supports three schemes for assigning the displacement threshold `Ed`:

* **`DTEMode = 1` (default) — uniform, species-wise.** Each element uses the
  single `Ed` given to its `Element(...)`.
* **`DTEMode = 2` — environment-dependent.** `Ed` is looked up from a table that
  maps the local defect environment (neighbouring vacancies) to an `Ed`,
  supplied via `DTEFile`. Useful for high-fluence 2D materials where pore-edge
  atoms have a different `Ed` from bulk atoms.
* **`DTEMode = 3` — fully custom.** `Ed` (and the binding energy) is returned by
  two user-defined functions, evaluated for every potential recoil:

  ```julia
  function GetDTECustom(atom::Atom, simulator::Simulator)
      # any rule of position, recoil direction, local environment, ...
      return Ed_value
  end
  GetBDECustom(atom::Atom, simulator::Simulator) = atom.bde
  ```

  Define them at top level in your script and set `DTEMode = 3`. The
  [`examples/1D_CNT_anisotropic-DTE`](examples/1D_CNT_anisotropic-DTE) case uses
  this to give a direction-dependent threshold (Ed = 14 eV for atoms knocked
  outward, 25 eV inward) on a carbon nanotube.

---

## 7. Examples

| Directory                                                                          | Demonstrates                                                                                   |
| ---------------------------------------------------------------------------------- | --------------------------------------------------------------------------------------------- |
| [`examples/3D_implantation_B-in-Si`](examples/3D_implantation_B-in-Si)             | 3D ion implantation, **dynamic loading**, depth profile `R_p`, tilt angle, thermal vibrations |
| [`examples/1D_CNT_anisotropic-DTE`](examples/1D_CNT_anisotropic-DTE)               | 1D nanotube, **custom direction-dependent `Ed`** (`DTEMode = 3`), sputtering-yield scan        |
| [`examples/load_structure_from_file_C60`](examples/load_structure_from_file_C60)   | loading an **arbitrary structure from a file** (static mode), per-collision dump              |

Older example scripts are kept under [`examples/archive/`](examples/archive/).

---

## 8. Output files

* `*.dump` — LAMMPS-style snapshots (id, type, x, y, z [, extra columns]),
  readable by OVITO/VMD.
* `*.csv` — tabulated observables written by `@record` (e.g. `R_p.csv`,
  `nV.csv`).

---

## 9. Citation

If you use ASIM, please cite the accompanying paper: *ASIM: An atomic-scale
simulation framework for ion irradiation on multi-dimensional materials*.
