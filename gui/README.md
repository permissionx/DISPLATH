# DISPLAΘ GUI

A simple graphical user interface for the DISPLAΘ ion irradiation simulator.

## Installation

First, install the required Julia packages:

```julia
using Pkg
Pkg.add(["Gtk", "JSON3"])
```

## Usage

### Running the GUI

```bash
cd DISPLATH/gui/
julia simple_gui.jl
```

Or from Julia REPL:
```julia
include("gui/simple_gui.jl")
main()
```

## Features

### Current Features

1. **Tabbed Interface**: Organized parameters into logical groups
   - **Material Tab**: Material selection, lattice parameters, box sizes
   - **Simulation Tab**: Physical parameters (pMax, temperature, etc.)
   - **Ion Tab**: Ion type, energy, and number of cascades

2. **Material Presets**: Built-in presets for common materials
   - SiC (Silicon Carbide)
   - hBN (hexagonal Boron Nitride) 
   - Graphene
   - Custom material option

3. **Ion Selection**: Support for different ion types
   - N (Nitrogen)
   - Ne (Neon)
   - Ar (Argon)
   - Kr (Krypton)
   - Xe (Xenon)

4. **Real-time Feedback**:
   - Progress bar during simulation
   - Status updates
   - Live output log
   - Results display

5. **Parameter Management**:
   - Save/load parameter sets as JSON files
   - Tooltips for parameter guidance
   - Automatic material defaults

### Key Parameters

- **Lattice Constant**: Material lattice parameter in Angstroms
- **Box Sizes**: Simulation cell dimensions (number of unit cells)
- **pMax**: Maximum impact parameter for collision detection
- **Vacancy Recover Distance**: Distance threshold for vacancy annihilation
- **Temperature**: Simulation temperature in Kelvin
- **Debye Temperature**: Material Debye temperature
- **Stop Energy**: Minimum energy threshold in eV
- **Ion Energy**: Incident ion energy in eV
- **Number of Cascades**: Total cascades to simulate

## Example Workflow

1. **Select Material**: Choose from SiC, hBN, Graphene, or Custom
2. **Set Parameters**: Adjust lattice constant, box sizes, and simulation parameters
3. **Choose Ion**: Select ion type and energy
4. **Run Simulation**: Click "Run Simulation" and monitor progress
5. **View Results**: Check output log for vacancy statistics
6. **Save Results**: Simulation automatically saves results to CSV file

## Output

The GUI generates:
- **Real-time log**: Progress updates and intermediate results
- **CSV file**: Final results with vacancy counts per cascade
- **Statistics**: Average vacancies, count ranges, completion status

## Extending the GUI

The GUI is designed to be easily extensible:

### Adding New Materials

1. Add material to `material_combo` in `create_gui()`
2. Update `update_material_defaults()` with material properties
3. Add material-specific setup in `run_simulation()`

### Adding New Ion Types

1. Add ion to `ion_combo` in `create_gui()`
2. Update `typeDict` in `run_simulation()` with ion properties

### Adding New Parameters

1. Add GUI widgets in appropriate tab
2. Update `get_parameters()` to include new parameters
3. Update save/load functions for persistence

## Current Limitations

1. **Material Support**: Currently only SiC is fully implemented
2. **File Dialogs**: Uses simple filename for save/load (no file browser)
3. **Visualization**: No real-time 3D visualization (planned for future)
4. **Advanced Features**: No DTE, KMC, or dynamic loading options yet

## Future Enhancements

- [ ] Full material library (hBN, Graphene, etc.)
- [ ] 3D visualization of simulation results
- [ ] Real-time plotting of vacancy statistics
- [ ] Advanced simulation modes (DTE, KMC)
- [ ] File browser dialogs
- [ ] Parameter validation
- [ ] Simulation pause/resume functionality
- [ ] Batch simulation capability

## Troubleshooting

### Common Issues

1. **Gtk not found**: Install with `Pkg.add("Gtk")`
2. **JSON3 not found**: Install with `Pkg.add("JSON3")`
3. **Simulation crashes**: Check parameter values and material setup
4. **GUI freezes**: Use smaller number of cascades for testing

### Performance Tips

1. Start with small box sizes (10x10x100) for testing
2. Use fewer cascades (10-50) during development
3. Monitor memory usage for large simulations
4. Save parameters before long runs