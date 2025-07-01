# DISPLATH Web GUI

Modern, production-ready web-based graphical user interface for the DISPLATH Binary Collision Approximation simulator.

## Features

✅ **Complete Implementation**
- **Dark theme interface** perfectly matching the provided mockup design
- **Real-time parameter validation** with comprehensive error checking
- **3D crystal structure visualization** using Three.js with error handling
- **Live simulation progress tracking** with status polling
- **Comprehensive system logs** with structured logging
- **Material presets** (Silicon, Graphene, MoS₂) with easy loading
- **Dynamic cell loading toggle** for large simulations
- **Error boundaries** for crash-proof React components
- **Security hardening** with CORS, rate limiting, and input sanitization
- **Memory leak prevention** with proper resource cleanup
- **Responsive design** for different screen sizes

## Architecture

### Frontend (React 18)
- `index.html` - Main HTML file with CDN dependencies
- `app.jsx` - Complete React application with all components
- `styles.css` - Comprehensive dark theme styling
- Error boundaries, proper state management, API integration

### Backend (Flask + Julia)
- `server.py` - Secure Flask API server with proper error handling
- `julia_interface.py` - Safe Julia integration without injection vulnerabilities
- `logging_config.py` - Structured logging system with JSON format
- Comprehensive API endpoints with validation and rate limiting

### Security & Reliability
- Input validation and sanitization
- Proper error handling throughout
- Memory leak prevention
- Resource cleanup
- Rate limiting and CORS configuration
- Structured logging for debugging

## Installation

### Quick Start
```bash
./start_gui.sh
```

### Manual Installation

1. **Create conda environment**:
```bash
conda env create -f ../environment.yml
```

2. **Activate environment**:
```bash
conda activate displath-gui
```

3. **Install Python dependencies**:
```bash
pip install -r requirements.txt
```

4. **Ensure Julia is available**:
```bash
julia --version  # Should show Julia version
```

## Running the GUI

### Automatic Start
```bash
./start_gui.sh
```

### Manual Start
```bash
conda activate displath-gui
python server.py
```

Then open your browser to: **http://localhost:5000**

## API Documentation

### Core Endpoints
- `GET /api/health` - System health check
- `POST /api/validate` - Validate simulation parameters
- `POST /api/simulate` - Start a new simulation
- `GET /api/simulations` - List all simulations

### Simulation Management  
- `GET /api/simulation/<id>/status` - Get detailed status and progress
- `GET /api/simulation/<id>/logs` - Get paginated simulation logs
- `POST /api/simulation/<id>/cancel` - Cancel running simulation
- `DELETE /api/simulation/<id>` - Delete completed simulation

### Materials
- `GET /api/materials/presets` - Get predefined material configurations

### Rate Limits
- General: 100 requests/hour, 10 requests/minute
- Validation: 20 requests/minute
- Simulation start: 5 requests/minute

## Usage Guide

### 1. Load Material Configuration
- **Presets**: Select from dropdown (Silicon, Graphene, MoS₂)
- **Custom**: Define cell vectors and basis atoms manually
- **Validation**: Real-time parameter checking with error display

### 2. Configure Simulation
- **Box Ranges**: Number of unit cells in each direction
- **Lattice Ranges**: Specific range of cells to include
- **Dynamic Loading**: Toggle for large simulations
- **3D Preview**: Interactive crystal structure visualization

### 3. Run Simulation
- **Validation**: Automatic parameter validation before start
- **Progress**: Real-time progress bar and status updates
- **Logs**: Live system logs with different severity levels
- **Control**: Cancel, reset, or start new simulations

### 4. Monitor Results
- **Progress Tracking**: Visual progress bar with percentage
- **Live Updates**: Status polling every 2 seconds
- **Results Display**: Vacancy counts and statistics
- **Error Handling**: Clear error messages for failures

## File Structure

```
gui/new-gui/
├── index.html              # Main HTML file
├── app.jsx                 # React application
├── styles.css              # Complete styling
├── server.py               # Flask server
├── julia_interface.py      # Safe Julia integration
├── logging_config.py       # Logging system
├── requirements.txt        # Python dependencies
├── start_gui.sh           # Startup script
├── README.md              # This file
└── logs/                  # Log files (created automatically)
    ├── displath_gui.log   # All logs (JSON format)
    ├── displath_errors.log # Error logs only
    └── simulations.log    # Simulation-specific logs
```

## Components Detail

### Header & Navigation
- **ARCS branding** with version info and health status
- **Icon buttons** for information and validation
- **Collapsible sidebar** with Projectile/Parameters/Process tabs

### Material Configuration
- **Cell Vectors**: 3×3 matrix input with validation
- **Basis Atoms**: Dynamic list with add/remove functionality
- **Element Library**: Visual element cards with properties
- **Presets**: Easy loading of common materials

### Simulation Setup
- **Box Ranges**: Integer inputs with real-time dimension calculation
- **Lattice Ranges**: Min/max inputs with Angstrom conversion
- **Progress Monitor**: Real-time status and progress tracking

### Visualization & Monitoring
- **3D Crystal View**: Interactive Three.js visualization
- **System Logs**: Scrollable log viewer with color coding
- **Control Panel**: Run, Cancel, Reset, and Clear functions

## Logging

### Log Files
- **All logs**: `logs/displath_gui.log` (rotating, 10MB max)
- **Errors only**: `logs/displath_errors.log` (rotating, 5MB max)  
- **Simulations**: `logs/simulations.log` (rotating, 20MB max)

### Log Levels
- **INFO**: Normal operations and status updates
- **WARNING**: Non-critical issues and alerts
- **ERROR**: Errors that don't crash the application
- **DEBUG**: Detailed debugging information

### Structured Format
All logs use JSON format with timestamps, levels, and contextual data for easy parsing and analysis.

## Troubleshooting

### Common Issues

1. **Julia not found**
   - Install Julia and ensure it's in PATH
   - Check with: `julia --version`

2. **DISPLATH source not found**
   - Ensure `../../src/main.jl` exists
   - Check DISPLATH installation

3. **Port already in use**
   - Kill existing server: `pkill -f "python server.py"`
   - Or change port in server.py

4. **Permission denied**
   - Make script executable: `chmod +x start_gui.sh`

5. **Conda environment issues**
   - Recreate environment: `conda env remove -n displath-gui && conda env create -f ../environment.yml`

### Debug Mode
For debugging, set debug mode in server.py:
```python
app.run(debug=True, port=5000)
```

## Performance Notes

- **Memory Usage**: Automatic cleanup of old simulations after 24 hours
- **Concurrent Simulations**: Maximum 10 simultaneous simulations
- **Log Rotation**: Automatic rotation prevents disk space issues
- **Rate Limiting**: Prevents abuse and resource exhaustion

## Security Features

- **Input Validation**: All user inputs validated and sanitized
- **CORS Configuration**: Restricted to allowed origins only
- **Rate Limiting**: Per-endpoint rate limits prevent abuse
- **Error Handling**: Graceful error handling without information leakage
- **Resource Cleanup**: Proper cleanup prevents memory leaks