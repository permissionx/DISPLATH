# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

DISPLATH (DISPLAΘ Is the Second Program Leveraging Accurate Θ) is a Julia-based Binary Collision Approximation (BCA) simulator for modeling ion irradiation processes in materials. It simulates collision cascades to predict radiation damage (vacancies) in materials.

## Key Commands

### Running Simulations
```bash
# Run a simulation example
julia examples/Static_load/graphene/main.jl

# Run with specific number of threads
julia -t 8 examples/Static_load/graphene/main.jl
```

### Environment Setup
```bash
# Run the installation script to set up ARCS environment
./install_arcs.sh

# The script sets up:
# - ARCS_REPO=$HOME/.arcs_repository
# - ARCS_HOME=/beegfs/science-share/arcs/DISPLATH
# - Julia PATH=/beegfs/science-share/julia/bin
```

### Testing
There is no formal test runner. Test files exist in the `test/` directory but are run individually:
```bash
julia test/test.jl
```

## Code Architecture

### Module Structure
The main module is loaded via:
```julia
include("/path/to/DISPLATH/src/DISPLATH.jl")
```

### Core Components
- **src/DISPLATH.jl**: Main module file that includes all components
- **src/bca.jl**: Binary Collision Approximation core logic
- **src/types.jl**: Type definitions (Atom, Simulator, Parameters)
- **src/geometry.jl**: Geometric calculations and grid operations
- **src/dynamics.jl**: Collision dynamics and cascade simulation
- **src/dte.jl**: Dynamic Time Evolution functionality
- **src/io.jl**: Input/output operations and data dumping

### Key Types
- `Atom`: Represents an atom with position, velocity, type
- `Simulator`: Main simulation container with atoms and grid
- `Parameters`: Simulation parameters including materials and thresholds

### Typical Simulation Workflow
1. Define material structure (lattice vectors, basis, atom types)
2. Create Parameters object with collision settings
3. Initialize Simulator with box size and grid
4. Create ions with specific energy/direction
5. Run Cascade! to simulate collisions
6. Analyze defects with DefectStatics

### Data Files
- **θτ files** (.thetatau): Scattering angle and time data
- **DTE files** (.dte): Dynamic time evolution parameters
- **Dump files** (.dump): Atom position snapshots

### Key Functions
- `Cascade!(atom, simulator)`: Run collision cascade
- `DefectStatics(simulator)`: Analyze vacancies and interstitials
- `Save!/Restore!`: State management for multiple runs
- `@dump`: Macro for saving atom positions
- `@record`: Macro for logging results to CSV

### Material Examples
The codebase includes examples for:
- 2D materials: graphene, hBN, MoS2
- 3D materials: Si, SiC
- Various ions: He, Ne, Ar, Kr, Xe

### Performance Considerations
- Uses Julia's threading for parallel execution
- StableRNG for reproducible parallel random numbers
- Static arrays for performance in geometric calculations

### Important Environment Variables
- `ARCS_HOME`: Points to DISPLATH installation directory
- `ARCS_REPO`: Repository for θτ and DTE data files
- `IS_DYNAMIC_LOAD`: Flag for dynamic vs static loading mode

## Development Notes

When modifying the simulator:
- Maintain thread safety when using parallel features
- Use `@dump` and `@record` macros for consistent output
- Follow existing patterns in examples/ for new simulations
- Defect analysis requires proper vacancy recovery distance settings