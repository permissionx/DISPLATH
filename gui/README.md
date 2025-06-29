# DISPLAΘ Web GUI

A modern web-based graphical user interface for the DISPLAΘ Binary Collision Approximation (BCA) simulator.

## Features

- **Modern Interface**: Clean, responsive design with tabbed navigation
- **Material Configuration**: Set up crystal structures, lattice parameters, and element types
- **Simulation Parameters**: Configure simulation settings including temperature, energy thresholds
- **Ion Configuration**: Define ion properties, energy, direction, and position
- **Static/Dynamic Loading**: Support for both static and dynamic cell loading modes
- **Real-time Results**: View simulation results with statistics and histogram plots
- **Preset Materials**: Pre-configured setups for common materials (Graphene, Silicon)

## Quick Start

1. **Start the Server**:
   ```julia
   cd gui
   julia server.jl
   ```

2. **Open Browser**: Navigate to `http://localhost:8080`

3. **Configure Simulation**:
   - Use preset materials or configure manually
   - Set simulation parameters
   - Define ion properties
   - Run simulation

## Interface Overview

### Material Tab
- **Preset Materials**: Quick setup for Graphene and Silicon
- **Primary Vectors**: Define unit cell lattice vectors
- **Lattice Ranges**: Set simulation box boundaries
- **Basis Atoms**: Configure atomic positions in unit cell
- **Element Types**: Define atomic species with displacement energies

### Simulation Tab
- **Box Sizes**: Simulation cell dimensions
- **Number of Runs**: Statistical sampling size
- **Physical Parameters**: Temperature, Debye temperature, pMax
- **Dynamic Loading**: Enable for large-scale simulations

### Ion Tab
- **Ion Properties**: Type, energy, direction, position
- **Energy Range**: 1 eV to 1 MeV supported
- **Direction Control**: 3D vector specification

### Output Tab
- **Results Display**: Statistics and histogram visualization
- **Data Export**: Raw simulation data access
- **Error Handling**: Clear error messages and debugging info

## Supported Materials

### Graphene
- Monolayer configuration
- C-C bond length: 1.42 Å
- Displacement energy: 19.96 eV

### Silicon
- Diamond cubic structure
- Lattice parameter: 5.431 Å
- Displacement energy: 20.0 eV (Si), 1.0 eV (B)

## Configuration Examples

### Basic Graphene Irradiation
```
Material: Graphene (preset)
Ion: Ne+ at 1000 eV
Direction: [0, 0, -1] (normal incidence)
Runs: 100
```

### Silicon Ion Implantation
```
Material: Silicon (preset)  
Ion: B+ at 10 keV
Dynamic Loading: Enabled
Runs: 1000
```

## Technical Details

### Backend Architecture
- **Julia HTTP Server**: Handles simulation requests
- **REST API**: JSON-based communication
- **CORS Support**: Cross-origin resource sharing enabled

### Frontend Stack
- **Vanilla JavaScript**: No framework dependencies
- **Chart.js**: Data visualization
- **Modern CSS**: Grid/Flexbox layouts
- **Responsive Design**: Mobile-friendly interface

### File Structure
```
gui/
├── server.jl          # Julia backend server
├── index.html         # Main GUI interface
├── style.css          # Modern styling
├── script.js          # Frontend logic
└── README.md          # This file
```

## API Endpoints

- `GET /`: Main interface
- `GET /api/presets`: Available material presets
- `POST /api/simulate`: Run simulation with parameters

## Troubleshooting

### Common Issues

1. **Server Won't Start**
   - Check Julia installation and dependencies
   - Verify port 8080 is available
   - Ensure DISPLAΘ core modules load correctly

2. **Simulation Errors**
   - Validate input parameters
   - Check thetatau files exist for element combinations
   - Verify material configuration is physically reasonable

3. **Performance Issues**
   - Reduce number of simulation runs
   - Enable dynamic loading for large systems
   - Check system memory availability

### Dependencies

Required Julia packages:
- HTTP.jl
- JSON3.jl
- DISPLAΘ core modules

## Development

### Adding New Materials
1. Create preset configuration in `server.jl`
2. Add to presets API endpoint
3. Update frontend preset selector

### Extending Functionality
- Modify `SimulationConfig` struct for new parameters
- Update form collection in `script.js`
- Add corresponding UI elements in `index.html`

## Performance Notes

- Recommended: 100-1000 runs for statistical significance
- Dynamic loading required for systems >10⁶ atoms
- Browser memory limits visualization to ~10⁴ data points

## Browser Compatibility

- Chrome/Chromium 90+
- Firefox 88+
- Safari 14+
- Edge 90+