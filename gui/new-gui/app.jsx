const { useState, useEffect, useRef, useCallback, useMemo } = React;

// Utility functions for safe type conversion
const safeParseFloat = (value, defaultValue = 0) => {
    const parsed = parseFloat(value);
    return isNaN(parsed) ? defaultValue : parsed;
};

const safeParseInt = (value, defaultValue = 0) => {
    const parsed = parseInt(value, 10);
    return isNaN(parsed) ? defaultValue : parsed;
};

// Error Boundary Component
class ErrorBoundary extends React.Component {
    constructor(props) {
        super(props);
        this.state = { hasError: false, error: null };
    }

    static getDerivedStateFromError(error) {
        return { hasError: true, error };
    }

    componentDidCatch(error, errorInfo) {
        console.error('GUI Error:', error, errorInfo);
    }

    render() {
        if (this.state.hasError) {
            return (
                <div className="error-boundary">
                    <h2>Something went wrong</h2>
                    <p>{this.state.error?.message}</p>
                    <button onClick={() => this.setState({ hasError: false, error: null })}>
                        Try again
                    </button>
                </div>
            );
        }

        return this.props.children;
    }
}

// API client with proper error handling
const apiClient = {
    baseURL: '',
    
    async request(endpoint, options = {}) {
        try {
            const response = await fetch(`${this.baseURL}/api${endpoint}`, {
                headers: {
                    'Content-Type': 'application/json',
                    ...options.headers
                },
                ...options
            });
            
            if (!response.ok) {
                const errorData = await response.json().catch(() => ({ error: 'Unknown error' }));
                throw new Error(errorData.error || `HTTP ${response.status}`);
            }
            
            return await response.json();
        } catch (error) {
            console.error(`API Error (${endpoint}):`, error);
            throw error;
        }
    },

    get(endpoint) {
        return this.request(endpoint, { method: 'GET' });
    },

    post(endpoint, data) {
        return this.request(endpoint, {
            method: 'POST',
            body: JSON.stringify(data)
        });
    },

    delete(endpoint) {
        return this.request(endpoint, { method: 'DELETE' });
    }
};

// 3D Crystal Visualization Component
const CrystalVisualization = ({ cellVectors, basisAtoms }) => {
    const containerRef = useRef(null);
    const sceneRef = useRef(null);
    const rendererRef = useRef(null);
    const animationRef = useRef(null);

    useEffect(() => {
        if (!containerRef.current || !window.THREE) return;

        // Clean up previous scene
        if (rendererRef.current) {
            rendererRef.current.dispose();
            containerRef.current.innerHTML = '';
        }

        try {
            // Scene setup
            const scene = new THREE.Scene();
            scene.background = new THREE.Color(0x333333);
            sceneRef.current = scene;

            // Camera
            const camera = new THREE.PerspectiveCamera(75, 300/250, 0.1, 1000);
            camera.position.set(3, 3, 3);
            camera.lookAt(0, 0, 0);

            // Renderer
            const renderer = new THREE.WebGLRenderer({ antialias: true });
            renderer.setSize(300, 250);
            renderer.shadowMap.enabled = true;
            renderer.shadowMap.type = THREE.PCFSoftShadowMap;
            rendererRef.current = renderer;

            containerRef.current.appendChild(renderer.domElement);

            // Lighting
            const ambientLight = new THREE.AmbientLight(0x404040, 0.4);
            scene.add(ambientLight);

            const directionalLight = new THREE.DirectionalLight(0xffffff, 0.8);
            directionalLight.position.set(10, 10, 5);
            directionalLight.castShadow = true;
            scene.add(directionalLight);

            // Create crystal structure
            const geometry = new THREE.BoxGeometry(2, 2, 0.1);
            const material = new THREE.MeshLambertMaterial({
                color: 0xd946ef,
                transparent: true,
                opacity: 0.8
            });

            const crystal = new THREE.Mesh(geometry, material);
            scene.add(crystal);

            // Add lattice lines
            const edgesGeometry = new THREE.EdgesGeometry(geometry);
            const edgesMaterial = new THREE.LineBasicMaterial({ color: 0xffffff, opacity: 0.3, transparent: true });
            const edges = new THREE.LineSegments(edgesGeometry, edgesMaterial);
            scene.add(edges);

            // Controls
            if (window.THREE.OrbitControls) {
                const controls = new THREE.OrbitControls(camera, renderer.domElement);
                controls.enableDamping = true;
                controls.dampingFactor = 0.05;
                controls.autoRotate = true;
                controls.autoRotateSpeed = 1.0;
            }

            // Animation loop
            const animate = () => {
                animationRef.current = requestAnimationFrame(animate);
                renderer.render(scene, camera);
            };
            animate();

        } catch (error) {
            console.error('Three.js error:', error);
            containerRef.current.innerHTML = `
                <div style="display: flex; flex-direction: column; align-items: center; justify-content: center; height: 100%; color: #999;">
                    <div>3D Visualization</div>
                    <div style="font-size: 11px; margin-top: 5px;">WebGL not available</div>
                </div>
            `;
        }

        return () => {
            if (animationRef.current) {
                cancelAnimationFrame(animationRef.current);
            }
            if (rendererRef.current) {
                rendererRef.current.dispose();
            }
            if (containerRef.current) {
                containerRef.current.innerHTML = '';
            }
        };
    }, [cellVectors, basisAtoms]);

    return <div ref={containerRef} style={{ width: '100%', height: '100%' }} />;
};

// Main Application Component
const DisplathGUI = () => {
    // State for sidebar and navigation
    const [sidebarCollapsed, setSidebarCollapsed] = useState(false);
    const [activeTab, setActiveTab] = useState('parameters');
    
    // Material and structure state
    const [cellVectors, setCellVectors] = useState([
        [4.260, 0.000, 0.000],
        [0.000, 4.263, 0.000], 
        [0.000, 0.000, 6.700]
    ]);
    
    const [basisAtoms, setBasisAtoms] = useState([
        { x: 0.0, y: 0.0, z: 0.0, type: 1 },
        { x: 0.0, y: 0.0, z: 0.0, type: 2 }
    ]);
    
    const [boxRanges, setBoxRanges] = useState([400, 400, 20]);
    
    const [latticeRanges, setLatticeRanges] = useState([
        { min: 0, max: 400 },
        { min: 0, max: 400 },
        { min: 2, max: 3 }
    ]);
    
    // UI state
    const [dynamicLoad, setDynamicLoad] = useState(true);
    const [logs, setLogs] = useState([
        { time: '00:37:03', message: 'System initialized. Ready for simulation.', level: 'info' },
        { time: '00:37:03', message: 'System initialized. Ready for simulation.', level: 'info' }
    ]);
    
    // Simulation state
    const [currentSimulation, setCurrentSimulation] = useState(null);
    const [validationErrors, setValidationErrors] = useState([]);
    const [materialPresets, setMaterialPresets] = useState({});

    // Computed values
    const dimensions = useMemo(() => ({
        x: cellVectors[0][0] * boxRanges[0],
        y: cellVectors[1][1] * boxRanges[1],
        z: cellVectors[2][2] * boxRanges[2]
    }), [cellVectors, boxRanges]);

    // Load material presets on component mount
    useEffect(() => {
        const loadPresets = async () => {
            try {
                const presets = await apiClient.get('/materials/presets');
                setMaterialPresets(presets);
            } catch (error) {
                console.error('Failed to load material presets:', error);
            }
        };
        
        loadPresets();
    }, []);

    // Material preset loader
    const loadMaterialPreset = useCallback((presetName) => {
        const preset = materialPresets[presetName];
        if (!preset) return;
        
        setCellVectors(preset.cellVectors);
        setBasisAtoms(preset.basisAtoms);
        setBoxRanges(preset.boxRanges);
        setLatticeRanges(preset.latticeRanges);
        
        addLog(`Loaded ${preset.name} material preset`, 'info');
    }, [materialPresets]);

    // Utility functions
    const addLog = useCallback((message, level = 'info') => {
        const time = new Date().toLocaleTimeString('en-US', { hour12: false });
        setLogs(prev => [...prev.slice(-49), { time, message, level }]);
    }, []);

    const clearLogs = useCallback(() => {
        setLogs([]);
    }, []);

    const resetSimulation = useCallback(() => {
        setCurrentSimulation(null);
        setValidationErrors([]);
        addLog('Simulation reset', 'info');
    }, [addLog]);

    const startSimulation = useCallback(async () => {
        try {
            addLog('Starting simulation...', 'info');
            
            const params = {
                cellVectors,
                basisAtoms,
                boxRanges,
                latticeRanges,
                dynamicLoad,
                numRuns: 100,
                ionEnergy: 100000,
                ionType: 1
            };

            // Validate parameters first
            const validation = await apiClient.post('/validate', params);
            
            if (!validation.valid) {
                setValidationErrors(validation.errors);
                addLog(`Validation failed: ${validation.errors.join(', ')}`, 'error');
                return;
            }

            setValidationErrors([]);
            
            // Start simulation
            const response = await apiClient.post('/simulate', params);
            setCurrentSimulation({ 
                id: response.simulationId, 
                status: 'running',
                progress: 0 
            });
            
            addLog(`Simulation started: ${response.simulationId}`, 'info');
            
        } catch (error) {
            addLog(`Failed to start simulation: ${error.message}`, 'error');
        }
    }, [cellVectors, basisAtoms, boxRanges, latticeRanges, dynamicLoad, addLog]);

    return (
        <ErrorBoundary>
            <div className="app">
                {/* Validation Errors */}
                {validationErrors.length > 0 && (
                    <div className="validation-errors">
                        <h4>Validation Errors:</h4>
                        <ul>
                            {validationErrors.map((error, index) => (
                                <li key={index}>{error}</li>
                            ))}
                        </ul>
                    </div>
                )}

                {/* Header */}
                <header className="header">
                    <div className="header-content">
                        <span className="logo">ARCS</span>
                        <span className="bullet">•</span>
                        <span className="subtitle">Atomic Resolution Collision Simulator Playground</span>
                        <span className="version">version 0.1</span>
                    </div>
                    <div className="header-actions">
                        <button className="icon-btn" title="Information">
                            <span>ℹ</span>
                        </button>
                        <button className="icon-btn" title="Export">
                            <span>🗃</span>
                        </button>
                    </div>
                </header>

                {/* Project Bar */}
                <div className="project-bar">
                    <div className="project-info">
                        <div className="project-title-row">
                            <span className="project-title">Graphene Sequential Irradiation</span>
                            <button className="btn btn-secondary">+ New Project</button>
                        </div>
                        <div className="project-path">
                            <span>📁</span>
                            <span>/home/Perry/researches/Graphene</span>
                        </div>
                    </div>
                    <div className="dynamic-load-toggle">
                        <span className="dynamic-load-label">Turn on Dynamic Load</span>
                        <label className="toggle-switch">
                            <input 
                                type="checkbox" 
                                checked={dynamicLoad} 
                                onChange={(e) => setDynamicLoad(e.target.checked)}
                            />
                            <span className="toggle-slider"></span>
                        </label>
                    </div>
                </div>

                {/* Main Layout */}
                <div className="main-layout">
                    {/* Sidebar */}
                    <div className={`sidebar ${sidebarCollapsed ? 'collapsed' : ''}`}>
                        <div className="sidebar-toggle">
                            <button 
                                className="icon-btn" 
                                onClick={() => setSidebarCollapsed(!sidebarCollapsed)}
                            >
                                <span>☰</span>
                            </button>
                        </div>
                        
                        <div 
                            className={`sidebar-item ${activeTab === 'projectile' ? 'active' : ''}`}
                            onClick={() => setActiveTab('projectile')}
                        >
                            <span className="sidebar-item-icon">🚀</span>
                            {!sidebarCollapsed && <span className="sidebar-item-text">Projectile</span>}
                        </div>
                        
                        <div 
                            className={`sidebar-item ${activeTab === 'parameters' ? 'active' : ''}`}
                            onClick={() => setActiveTab('parameters')}
                        >
                            <span className="sidebar-item-icon">🔧</span>
                            {!sidebarCollapsed && <span className="sidebar-item-text">Parameters</span>}
                        </div>
                        
                        <div 
                            className={`sidebar-item ${activeTab === 'process' ? 'active' : ''}`}
                            onClick={() => setActiveTab('process')}
                        >
                            <span className="sidebar-item-icon">⚡</span>
                            {!sidebarCollapsed && <span className="sidebar-item-text">Process</span>}
                        </div>
                    </div>

                    {/* Content Area */}
                    <div className="content-area">
                        {/* Left Column - Parameters */}
                        <div className="column left-column">
                            {/* Cell Vectors Section */}
                            <div className="section">
                                <div className="section-header">
                                    <h3 className="section-title">
                                        Cell Vectors (Å)
                                        <span className="info-icon" title="Unit cell vectors in Angstroms">ℹ</span>
                                    </h3>
                                </div>
                                <div className="cell-vectors-matrix">
                                    {cellVectors.map((row, i) => (
                                        <div key={`row-${i}`} className="matrix-row">
                                            {row.map((value, j) => (
                                                <input
                                                    key={`cv-${i}-${j}`}
                                                    type="number"
                                                    className="input-field matrix-input"
                                                    value={value}
                                                    step="0.001"
                                                    onChange={(e) => {
                                                        const newVectors = [...cellVectors];
                                                        newVectors[i][j] = safeParseFloat(e.target.value);
                                                        setCellVectors(newVectors);
                                                    }}
                                                />
                                            ))}
                                        </div>
                                    ))}
                                </div>
                            </div>

                            {/* Basis Atoms Section */}
                            <div className="section">
                                <div className="section-header">
                                    <h3 className="section-title">Basis Atoms</h3>
                                </div>
                                <div className="basis-atoms-container">
                                    {basisAtoms.map((atom, index) => (
                                        <div key={index} className="basis-atom-row">
                                            <input
                                                type="number"
                                                className="input-field basis-input"
                                                value={atom.x}
                                                step="0.001"
                                                onChange={(e) => {
                                                    const newAtoms = [...basisAtoms];
                                                    newAtoms[index].x = safeParseFloat(e.target.value);
                                                    setBasisAtoms(newAtoms);
                                                }}
                                            />
                                            <input
                                                type="number"
                                                className="input-field basis-input"
                                                value={atom.y}
                                                step="0.001"
                                                onChange={(e) => {
                                                    const newAtoms = [...basisAtoms];
                                                    newAtoms[index].y = safeParseFloat(e.target.value);
                                                    setBasisAtoms(newAtoms);
                                                }}
                                            />
                                            <input
                                                type="number"
                                                className="input-field basis-input"
                                                value={atom.z}
                                                step="0.001"
                                                onChange={(e) => {
                                                    const newAtoms = [...basisAtoms];
                                                    newAtoms[index].z = safeParseFloat(e.target.value);
                                                    setBasisAtoms(newAtoms);
                                                }}
                                            />
                                            <input
                                                type="number"
                                                className="input-field basis-type"
                                                value={atom.type}
                                                min="1"
                                                onChange={(e) => {
                                                    const newAtoms = [...basisAtoms];
                                                    newAtoms[index].type = safeParseInt(e.target.value, 1);
                                                    setBasisAtoms(newAtoms);
                                                }}
                                            />
                                            <button
                                                className="delete-btn"
                                                onClick={() => {
                                                    const newAtoms = basisAtoms.filter((_, i) => i !== index);
                                                    setBasisAtoms(newAtoms);
                                                }}
                                            >
                                                ×
                                            </button>
                                        </div>
                                    ))}
                                </div>
                                <button 
                                    className="btn btn-add"
                                    onClick={() => setBasisAtoms([...basisAtoms, { x: 0, y: 0, z: 0, type: 1 }])}
                                >
                                    + New Basis Atom
                                </button>
                            </div>

                            {/* Element Library Section */}
                            <div className="section">
                                <div className="section-header">
                                    <h3 className="section-title">Element Library</h3>
                                </div>
                                <div className="element-library">
                                    <div className="element-row">
                                        <div className="element-card">
                                            <div className="element-number">1</div>
                                            <div className="element-symbol">Si</div>
                                            <div className="element-properties">
                                                <div>DTE 22.2 eV</div>
                                                <div>BDE 11.1 eV</div>
                                                <div>Z    23</div>
                                                <div>Mass 23.2</div>
                                                <div>R    23.2</div>
                                            </div>
                                        </div>
                                        
                                        <div className="element-card">
                                            <div className="element-number">2</div>
                                            <div className="element-symbol">C</div>
                                            <div className="element-properties">
                                                <div>DTE 22.2 eV</div>
                                                <div>BDE 11.1 eV</div>
                                                <div>Z    23</div>
                                                <div>Mass 23.2</div>
                                                <div>R    23.2</div>
                                            </div>
                                        </div>
                                        
                                        <div className="element-card">
                                            <div className="element-number">3</div>
                                            <div className="element-symbol">N</div>
                                            <div className="element-properties">
                                                <div>DTE 22.2 eV</div>
                                                <div>BDE 11.1 eV</div>
                                                <div>Z    23</div>
                                                <div>Mass 23.2</div>
                                                <div>R    23.2</div>
                                            </div>
                                        </div>
                                    </div>
                                </div>
                                <button className="btn btn-add">
                                    + New Element
                                </button>
                            </div>
                        </div>

                        {/* Middle Column - Ranges */}
                        <div className="column middle-column">
                            {/* Box Ranges Section */}
                            <div className="section">
                                <div className="section-header">
                                    <h3 className="section-title">
                                        Box Ranges
                                        <span className="info-icon" title="Number of unit cells in each direction">ℹ</span>
                                    </h3>
                                </div>
                                <div className="box-ranges-inputs">
                                    {boxRanges.map((value, index) => (
                                        <input
                                            key={index}
                                            type="number"
                                            className="input-field box-input"
                                            value={value}
                                            min="1"
                                            onChange={(e) => {
                                                const newRanges = [...boxRanges];
                                                newRanges[index] = safeParseInt(e.target.value, 1);
                                                setBoxRanges(newRanges);
                                            }}
                                        />
                                    ))}
                                </div>
                                <div className="dimensions-info">
                                    <div>Dimensions: {dimensions.x.toFixed(1)} Å × {dimensions.y.toFixed(1)} Å × {dimensions.z.toFixed(1)} Å</div>
                                    <div>Center Point: ({(dimensions.x/2).toFixed(1)}, {(dimensions.y/2).toFixed(1)}, {(dimensions.z/2).toFixed(1)})</div>
                                </div>
                            </div>

                            {/* Lattice Range Section */}
                            <div className="section">
                                <div className="section-header">
                                    <h3 className="section-title">
                                        Lattice Range
                                        <span className="info-icon" title="Range of cells to include in simulation">ℹ</span>
                                    </h3>
                                </div>
                                <div className="lattice-ranges-container">
                                    {latticeRanges.map((range, index) => (
                                        <div key={index} className="lattice-range-row">
                                            <input
                                                type="number"
                                                className="input-field lattice-input"
                                                value={range.min}
                                                onChange={(e) => {
                                                    const newRanges = [...latticeRanges];
                                                    newRanges[index].min = safeParseInt(e.target.value, 0);
                                                    setLatticeRanges(newRanges);
                                                }}
                                            />
                                            <span className="lattice-to">to</span>
                                            <input
                                                type="number"
                                                className="input-field lattice-input"
                                                value={range.max}
                                                onChange={(e) => {
                                                    const newRanges = [...latticeRanges];
                                                    newRanges[index].max = safeParseInt(e.target.value, 1);
                                                    setLatticeRanges(newRanges);
                                                }}
                                            />
                                            <div className="lattice-range-value">
                                                {(range.min * cellVectors[index][index]).toFixed(1)} Å
                                            </div>
                                            <div className="lattice-range-value">
                                                {(range.max * cellVectors[index][index]).toFixed(1)} Å
                                            </div>
                                        </div>
                                    ))}
                                </div>
                            </div>
                        </div>

                        {/* Right Column - Visualization */}
                        <div className="column right-column">
                            {/* Crystal Visualization */}
                            <div className="section">
                                <div className="section-header">
                                    <h3 className="section-title">
                                        Crystal Visualization
                                        <span className="info-icon" title="3D view of crystal structure">ℹ</span>
                                    </h3>
                                </div>
                                <div className="crystal-viz">
                                    <CrystalVisualization 
                                        cellVectors={cellVectors} 
                                        basisAtoms={basisAtoms} 
                                    />
                                </div>
                                <div className="crystal-controls">
                                    <button className="btn btn-secondary">Reset</button>
                                </div>
                            </div>

                            {/* System Logs */}
                            <div className="section">
                                <div className="section-header">
                                    <h3 className="section-title">System Logs</h3>
                                </div>
                                <div className="system-logs">
                                    {logs.map((log, index) => (
                                        <div key={index} className={`log-entry log-${log.level}`}>
                                            <span className="log-time">[{log.time}]</span> {log.message}
                                        </div>
                                    ))}
                                </div>
                            </div>
                            
                            {/* Control Buttons */}
                            <div className="control-buttons">
                                <button className="btn btn-secondary" onClick={clearLogs}>
                                    Clear
                                </button>
                                <button className="btn btn-secondary" onClick={resetSimulation}>
                                    🔄 Reset
                                </button>
                                <button 
                                    className={`btn btn-primary ${currentSimulation?.status === 'running' ? 'loading' : ''}`}
                                    onClick={startSimulation}
                                    disabled={currentSimulation?.status === 'running'}
                                >
                                    ▶ Run Simulation
                                </button>
                            </div>
                        </div>
                    </div>
                </div>
            </div>
        </ErrorBoundary>
    );
};

// Render the application
const container = document.getElementById('app');
const root = ReactDOM.createRoot(container);
root.render(<DisplathGUI />);