# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Overview

DISPLAΘ is a Julia-based Binary Collision Approximation (BCA) simulator for modeling ion irradiation processes in materials. It simulates collision cascades to predict damage like vacancies produced by ions at various energies.

## Core Architecture

### Main Components

- **`src/main.jl`**: Entry point that includes all modules and dependencies
- **`src/types.jl`**: Core data structures (Atom, Simulator, Parameters, GridCell, etc.)
- **`src/bca.jl`**: Binary Collision Approximation physics calculations 
- **`src/dynamics.jl`**: Target finding and collision dynamics
- **`src/geometry.jl`**: Spatial geometry operations and atom creation
- **`src/dynamic_load.jl`**: Dynamic cell loading for large simulations
- **`src/dte.jl`**: Dynamic Time Evolution simulation capabilities
- **`src/kmc.jl`**: Kinetic Monte Carlo methods
- **`src/io.jl`**: File I/O and data output
- **`src/utils.jl`**: Utility functions

### Key Data Structures

- **`Simulator`**: Main simulation container holding atoms, lattice points, cell grid, and parameters
- **`Atom`**: Represents individual atoms with position, velocity, energy, type, and collision state
- **`GridCell`**: Spatial partitioning cells containing atoms and lattice points
- **`Parameters`**: Configuration container for simulation settings, material properties, and physical constants

### Simulation Flow

1. **Initialization**: Create simulator with lattice structure and material parameters
2. **Ion Injection**: Add ion with specified energy and direction  
3. **Cascade Simulation**: Use BCA to simulate collision cascades through `Cascade!()` function
4. **Target Finding**: Use spatial grid to efficiently find collision targets
5. **Analysis**: Count vacancies and output results

## Common Development Tasks

### Running Simulations

All simulations are executed by running Julia scripts that include `src/main.jl`:

```julia
include("src/main.jl")
# Define parameters and run simulation
```

Examples are in the `examples/` directory with working scripts like:
- `examples/SiC/run.jl`: SiC irradiation simulation
- `examples/B-Si/run.jl`: B on Si simulation

### Performance Analysis

Use the performance analysis tools in `examples/SiC/`:
- `quick_profile.jl`: Fast performance check
- `loadcell_profile.jl`: Analyze LoadCell function bottlenecks  
- `profile_analysis.jl`: Comprehensive performance analysis

### Material Setup

Material properties are defined through:
- `primaryVectors`: Unit cell vectors
- `basis`: Atomic positions in unit cell  
- `basisTypes`: Atom types for each basis position
- `typeDict`: Element properties (mass, radius, displacement energies)
- Theta-tau files in `thetatau_repository/` for scattering calculations

### Key Functions

- `Cascade!(ion, simulator)`: Main cascade simulation
- `LoadCell!()`: Dynamic loading of simulation cells
- `ShotTarget()`: Find collision targets using spatial grid
- `CountVacancies()`: Count defects after simulation
- `Save!()/Restore!()`: State management for multiple runs

### File Formats

- `.dump`: Atomic structure files (similar to LAMMPS format)
- `.dte`: Dynamic Time Evolution data files
- `.thetatau`: Scattering angle and energy transfer data
- `.csv`: Results output (vacancy counts, statistics)

## Dependencies

Required Julia packages:
- LinearAlgebra
- QuadGK  
- ProgressMeter
- StableRNGs
- Interpolations
- Dates

Optional for analysis:
- BenchmarkTools (for performance profiling)

## Directory Structure

- `src/`: Core simulation code
- `examples/`: Working simulation examples organized by material system
- `test/`: Test scripts and validation
- `thetatau_repository/`: Scattering data files
- `dte_repository/`: Dynamic Time Evolution data
- `dumps/`: Example atomic structure files
- `my_calculations/`: Research calculations and benchmarks

## Performance Considerations

- The `LoadCell!()` function for dynamic loading is often a bottleneck
- Spatial grid (`GridCell`) system enables efficient neighbor finding
- Memory management is critical for large simulations with many atoms
- Multi-threading support available through `Base.Threads`