const { useState, useEffect, useRef, useCallback, useMemo } = React;

// Utility functions for safe type conversion
const safeParseFloat = (value, defaultValue = 0) => {
    const parsed = parseFloat(value);
    return isNaN(parsed) ? defaultValue : parsed;
};

const safeParseInt = (value, defaultValue = 0) => {
    const parsed = parseInt(value);
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
        console.error('React Error Boundary caught an error:', error, errorInfo);
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

// Toggle Switch Component
const ToggleSwitch = ({ checked, onChange, disabled = false }) => {
    return (
        <label className="toggle-switch">
            <input 
                type="checkbox" 
                checked={checked} 
                onChange={onChange} 
                disabled={disabled}
            />
            <span className="toggle-slider"></span>
        </label>
    );
};

// Element Card Component
const ElementCard = ({ number, symbol, properties, onClick }) => {
    return (
        <div className="element-card" onClick={onClick} style={{ cursor: onClick ? 'pointer' : 'default' }}>
            <div className="element-number">{number}</div>
            <div className="element-symbol">{symbol}</div>
            <div className="element-properties">
                {Object.entries(properties).map(([key, value]) => (
                    <div key={key}>{key} {value}</div>
                ))}
            </div>
        </div>
    );
};

// Basis Atom Row Component with proper validation
const BasisAtomRow = ({ atom, onUpdate, onDelete, index }) => {
    const handleUpdate = useCallback((field, value) => {
        const numericValue = field === 'type' ? safeParseInt(value, 1) : safeParseFloat(value, 0);
        onUpdate({ ...atom, [field]: numericValue });
    }, [atom, onUpdate]);

    return (
        <div className="basis-atom-row">
            <input
                type="number"
                className="input-field"
                value={atom.x}
                onChange={(e) => handleUpdate('x', e.target.value)}
                step="0.001"
                placeholder="X"
            />
            <input
                type="number"
                className="input-field"
                value={atom.y}
                onChange={(e) => handleUpdate('y', e.target.value)}
                step="0.001"
                placeholder="Y"
            />
            <input
                type="number"
                className="input-field"
                value={atom.z}
                onChange={(e) => handleUpdate('z', e.target.value)}
                step="0.001"
                placeholder="Z"
            />
            <input
                type="number"
                className="input-field basis-atom-type"
                value={atom.type}
                onChange={(e) => handleUpdate('type', e.target.value)}
                min="1"
                placeholder="Type"
            />
            <button className="delete-btn" onClick={onDelete} title={`Delete atom ${index + 1}`}>
                ×
            </button>
        </div>
    );
};

// API client for backend communication
const API_BASE = '';

const apiClient = {
    async request(endpoint, options = {}) {
        try {
            const response = await fetch(`${API_BASE}${endpoint}`, {
                headers: {
                    'Content-Type': 'application/json',
                    ...options.headers
                },
                ...options
            });

            if (!response.ok) {
                const errorData = await response.json().catch(() => ({}));
                throw new Error(errorData.error || `HTTP ${response.status}`);
            }

            return await response.json();
        } catch (error) {
            console.error(`API request failed: ${endpoint}`, error);
            throw error;
        }
    },

    async validate(params) {
        return this.request('/api/validate', {
            method: 'POST',
            body: JSON.stringify(params)
        });
    },

    async startSimulation(params) {
        return this.request('/api/simulate', {
            method: 'POST',
            body: JSON.stringify(params)
        });
    },

    async getStatus(simulationId) {
        return this.request(`/api/simulation/${simulationId}/status`);
    },

    async cancelSimulation(simulationId) {
        return this.request(`/api/simulation/${simulationId}/cancel`, {
            method: 'POST'
        });
    },

    async getPresets() {
        return this.request('/api/materials/presets');
    },

    async healthCheck() {
        return this.request('/api/health');
    }
};

// Crystal Visualization Component with error handling
const CrystalVisualization = ({ cellVectors, basisAtoms }) => {
    const mountRef = useRef(null);
    const sceneRef = useRef(null);
    const rendererRef = useRef(null);
    const cameraRef = useRef(null);
    const animationIdRef = useRef(null);
    const [error, setError] = useState(null);

    useEffect(() => {
        if (!mountRef.current) return;
        
        try {
            // Check WebGL support
            if (!window.WebGLRenderingContext) {
                throw new Error('WebGL not supported');
            }

            // Scene setup
            const scene = new THREE.Scene();
            scene.background = new THREE.Color(0x3a3a3a);
            sceneRef.current = scene;

            // Camera setup
            const camera = new THREE.PerspectiveCamera(
                45,
                mountRef.current.clientWidth / mountRef.current.clientHeight,
                0.1,
                1000
            );
            camera.position.set(5, 5, 5);
            camera.lookAt(0, 0, 0);
            cameraRef.current = camera;

            // Renderer setup with error handling
            const renderer = new THREE.WebGLRenderer({ 
                antialias: true,
                alpha: true,
                preserveDrawingBuffer: true
            });
            
            renderer.setSize(mountRef.current.clientWidth, mountRef.current.clientHeight);
            renderer.shadowMap.enabled = true;
            renderer.shadowMap.type = THREE.PCFSoftShadowMap;
            rendererRef.current = renderer;
            
            mountRef.current.appendChild(renderer.domElement);

            // Create crystal structure visualization
            updateCrystalStructure(scene, cellVectors, basisAtoms);

            // Add lighting
            const ambientLight = new THREE.AmbientLight(0x404040, 0.6);
            scene.add(ambientLight);
            
            const directionalLight = new THREE.DirectionalLight(0xffffff, 0.8);
            directionalLight.position.set(5, 5, 5);
            directionalLight.castShadow = true;
            scene.add(directionalLight);

            // Add grid
            const gridHelper = new THREE.GridHelper(10, 10, 0x666666, 0x333333);
            scene.add(gridHelper);

            // Animation
            const animate = () => {
                if (!mountRef.current) return;
                
                animationIdRef.current = requestAnimationFrame(animate);
                
                // Slow rotation
                scene.rotation.y += 0.002;
                
                try {
                    renderer.render(scene, camera);
                } catch (renderError) {
                    console.error('Render error:', renderError);
                    setError('Rendering failed');
                    return;
                }
            };
            animate();

            // Handle resize
            const handleResize = () => {
                if (!mountRef.current || !camera || !renderer) return;
                
                camera.aspect = mountRef.current.clientWidth / mountRef.current.clientHeight;
                camera.updateProjectionMatrix();
                renderer.setSize(mountRef.current.clientWidth, mountRef.current.clientHeight);
            };
            
            window.addEventListener('resize', handleResize);

            return () => {
                window.removeEventListener('resize', handleResize);
                
                if (animationIdRef.current) {
                    cancelAnimationFrame(animationIdRef.current);
                }
                
                if (scene) {
                    scene.clear();
                }
                
                if (renderer) {
                    renderer.dispose();
                }
                
                if (mountRef.current && renderer?.domElement) {
                    try {
                        mountRef.current.removeChild(renderer.domElement);
                    } catch (e) {
                        // Element might already be removed
                    }
                }
            };

        } catch (error) {
            console.error('Three.js initialization error:', error);
            setError(error.message);
        }
    }, [cellVectors, basisAtoms]);

    const updateCrystalStructure = (scene, cellVectors, basisAtoms) => {
        // Remove existing atoms
        const existingAtoms = scene.children.filter(child => child.userData?.isAtom);
        existingAtoms.forEach(atom => scene.remove(atom));

        if (!cellVectors || !basisAtoms) return;

        // Create unit cell outline
        const cellGeometry = new THREE.BoxGeometry(
            cellVectors[0][0] || 1,
            cellVectors[1][1] || 1, 
            cellVectors[2][2] || 1
        );
        const cellMaterial = new THREE.MeshBasicMaterial({ 
            color: 0x4CAF50,
            wireframe: true,
            opacity: 0.6,
            transparent: true
        });
        const cellMesh = new THREE.Mesh(cellGeometry, cellMaterial);
        cellMesh.userData.isAtom = false;
        scene.add(cellMesh);

        // Add basis atoms
        basisAtoms.forEach((atom, index) => {
            const atomGeometry = new THREE.SphereGeometry(0.3, 16, 12);
            const atomMaterial = new THREE.MeshPhongMaterial({ 
                color: getAtomColor(atom.type)
            });
            const atomMesh = new THREE.Mesh(atomGeometry, atomMaterial);
            
            atomMesh.position.set(
                safeParseFloat(atom.x) * (cellVectors[0][0] || 1),
                safeParseFloat(atom.y) * (cellVectors[1][1] || 1),
                safeParseFloat(atom.z) * (cellVectors[2][2] || 1)
            );
            
            atomMesh.userData.isAtom = true;
            atomMesh.userData.atomType = atom.type;
            scene.add(atomMesh);
        });
    };

    const getAtomColor = (type) => {
        const colors = {
            1: 0x4CAF50,  // Si - Green
            2: 0x2196F3,  // C - Blue  
            3: 0xFF9800,  // N - Orange
        };
        return colors[type] || 0xFFFFFF;
    };

    if (error) {
        return (
            <div className="crystal-viz-error">
                <p>Visualization Error: {error}</p>
                <button onClick={() => setError(null)}>Retry</button>
            </div>
        );
    }

    return <div id="crystal-container" ref={mountRef}></div>;
};

// Main App Component with complete API integration
const App = () => {
    // UI State
    const [sidebarCollapsed, setSidebarCollapsed] = useState(false);
    const [dynamicLoad, setDynamicLoad] = useState(true);
    const [activeTab, setActiveTab] = useState('projectile');
    
    // Simulation State
    const [simulationId, setSimulationId] = useState(null);
    const [simulationStatus, setSimulationStatus] = useState(null);
    const [isSimulating, setIsSimulating] = useState(false);
    const [validationErrors, setValidationErrors] = useState([]);
    
    // Data State with proper typing
    const [cellVectors, setCellVectors] = useState([
        [4.260, 0.000, 0.000],
        [0.000, 4.263, 0.000],
        [0.000, 0.000, 6.700]
    ]);
    
    const [boxRanges, setBoxRanges] = useState([400, 400, 20]);
    
    const [latticeRanges, setLatticeRanges] = useState([
        { min: 0, max: 400 },
        { min: 0, max: 400 },
        { min: 2, max: 3 }
    ]);
    
    const [basisAtoms, setBasisAtoms] = useState([
        { id: generateId(), x: 0.0, y: 0.0, z: 0.0, type: 1 },
        { id: generateId(), x: 0.333, y: 0.667, z: 0.0, type: 2 }
    ]);
    
    const [elements] = useState([
        { number: 1, symbol: 'Si', properties: { DTE: '22.2 eV', BDE: '11.1 eV', Z: '14', Mass: '28.09', R: '1.46' } },
        { number: 2, symbol: 'C', properties: { DTE: '22.2 eV', BDE: '11.1 eV', Z: '6', Mass: '12.01', R: '1.46' } },
        { number: 3, symbol: 'N', properties: { DTE: '22.2 eV', BDE: '11.1 eV', Z: '7', Mass: '14.01', R: '1.46' } }
    ]);
    
    const [logs, setLogs] = useState([
        { time: 'Ready', message: 'System initialized. Ready for simulation.', level: 'info' },
        { time: new Date().toLocaleTimeString(), message: 'GUI loaded successfully.', level: 'info' }
    ]);
    
    const [materialPresets, setMaterialPresets] = useState({});
    const [healthStatus, setHealthStatus] = useState(null);

    // Helper function to generate unique IDs
    function generateId() {
        return Math.random().toString(36).substr(2, 9);
    }

    // Calculate dimensions with proper error handling
    const dimensions = useMemo(() => {
        try {
            return {
                x: (cellVectors[0][0] * boxRanges[0]).toFixed(1),
                y: (cellVectors[1][1] * boxRanges[1]).toFixed(1),
                z: (cellVectors[2][2] * boxRanges[2]).toFixed(1)
            };
        } catch (error) {
            console.error('Error calculating dimensions:', error);
            return { x: '0.0', y: '0.0', z: '0.0' };
        }
    }, [cellVectors, boxRanges]);
    
    const centerPoint = useMemo(() => {
        try {
            return {
                x: (parseFloat(dimensions.x) / 2).toFixed(1),
                y: (parseFloat(dimensions.y) / 2).toFixed(1),
                z: (parseFloat(dimensions.z) / 2).toFixed(1)
            };
        } catch (error) {
            console.error('Error calculating center point:', error);
            return { x: '0.0', y: '0.0', z: '0.0' };
        }
    }, [dimensions]);

    // Load presets and health check on mount
    useEffect(() => {
        loadMaterialPresets();
        checkHealth();
        
        // Set up periodic health check
        const healthInterval = setInterval(checkHealth, 30000); // Every 30 seconds
        return () => clearInterval(healthInterval);
    }, []);

    // Poll simulation status
    useEffect(() => {
        if (simulationId && isSimulating) {
            const pollInterval = setInterval(async () => {
                try {
                    const status = await apiClient.getStatus(simulationId);
                    setSimulationStatus(status);
                    
                    // Update logs if available
                    if (status.logs) {
                        setLogs(prevLogs => [...prevLogs, ...status.logs]);
                    }
                    
                    // Stop polling if simulation finished
                    if (status.status === 'completed' || status.status === 'error' || status.status === 'cancelled') {
                        setIsSimulating(false);
                        addLog(`Simulation ${status.status}`, status.status === 'completed' ? 'info' : 'error');
                    }
                } catch (error) {
                    console.error('Failed to get simulation status:', error);
                    addLog(`Failed to get simulation status: ${error.message}`, 'error');
                }
            }, 2000); // Poll every 2 seconds
            
            return () => clearInterval(pollInterval);
        }
    }, [simulationId, isSimulating]);

    // Helper functions
    const addLog = useCallback((message, level = 'info') => {
        const timestamp = new Date().toLocaleTimeString();
        setLogs(prevLogs => [...prevLogs, { time: timestamp, message, level }]);
    }, []);

    const loadMaterialPresets = async () => {
        try {
            const presets = await apiClient.getPresets();
            setMaterialPresets(presets);
            addLog('Material presets loaded');
        } catch (error) {
            console.error('Failed to load material presets:', error);
            addLog(`Failed to load presets: ${error.message}`, 'error');
        }
    };

    const checkHealth = async () => {
        try {
            const health = await apiClient.healthCheck();
            setHealthStatus(health);
            
            if (!health.julia_available) {
                addLog('Warning: Julia not available - simulations will fail', 'warning');
            }
        } catch (error) {
            console.error('Health check failed:', error);
            setHealthStatus({ status: 'unhealthy', julia_available: false });
            addLog('Server health check failed', 'error');
        }
    };

    const buildSimulationParams = () => {
        return {
            cellVectors,
            boxRanges,
            latticeRanges,
            basisAtoms: basisAtoms.map(atom => ({
                x: safeParseFloat(atom.x),
                y: safeParseFloat(atom.y),
                z: safeParseFloat(atom.z),
                type: safeParseInt(atom.type, 1)
            })),
            dynamicLoad,
            numRuns: 10, // Default small number for testing
            ionType: 1,
            ionEnergy: 100000.0,
            temperature: 300.0,
            stopEnergy: 20.0
        };
    };

    // Event handlers with proper validation
    const updateCellVector = useCallback((i, j, value) => {
        const numValue = safeParseFloat(value, 0);
        setCellVectors(prev => {
            const newVectors = prev.map(row => [...row]);
            newVectors[i][j] = numValue;
            return newVectors;
        });
    }, []);
    
    const updateBoxRange = useCallback((index, value) => {
        const numValue = safeParseInt(value, 1);
        setBoxRanges(prev => {
            const newRanges = [...prev];
            newRanges[index] = numValue;
            return newRanges;
        });
    }, []);
    
    const updateLatticeRange = useCallback((index, field, value) => {
        const numValue = safeParseInt(value, 0);
        setLatticeRanges(prev => {
            const newRanges = [...prev];
            newRanges[index] = { ...newRanges[index], [field]: numValue };
            return newRanges;
        });
    }, []);
    
    const addBasisAtom = useCallback(() => {
        const newAtom = {
            id: generateId(),
            x: 0.0,
            y: 0.0,
            z: 0.0,
            type: 1
        };
        setBasisAtoms(prev => [...prev, newAtom]);
        addLog('New basis atom added');
    }, [addLog]);
    
    const updateBasisAtom = useCallback((id, updatedAtom) => {
        setBasisAtoms(prev => prev.map(atom => 
            atom.id === id ? { ...updatedAtom, id } : atom
        ));
    }, []);
    
    const deleteBasisAtom = useCallback((id) => {
        setBasisAtoms(prev => prev.filter(atom => atom.id !== id));
        addLog('Basis atom removed');
    }, [addLog]);

    const validateParameters = async () => {
        try {
            const params = buildSimulationParams();
            const validation = await apiClient.validate(params);
            
            setValidationErrors(validation.errors || []);
            
            if (validation.valid) {
                addLog('Parameters validated successfully');
            } else {
                addLog(`Validation failed: ${validation.errors.join(', ')}`, 'error');
            }
            
            return validation.valid;
        } catch (error) {
            console.error('Validation failed:', error);
            addLog(`Validation error: ${error.message}`, 'error');
            return false;
        }
    };
    
    const runSimulation = async () => {
        if (isSimulating) {
            addLog('Simulation already running', 'warning');
            return;
        }
        
        addLog('Validating parameters...');
        const isValid = await validateParameters();
        
        if (!isValid) {
            addLog('Parameter validation failed', 'error');
            return;
        }
        
        try {
            setIsSimulating(true);
            addLog('Starting simulation...');
            
            const params = buildSimulationParams();
            const response = await apiClient.startSimulation(params);
            
            setSimulationId(response.simulationId);
            addLog(`Simulation started with ID: ${response.simulationId}`);
            
        } catch (error) {
            console.error('Failed to start simulation:', error);
            addLog(`Failed to start simulation: ${error.message}`, 'error');
            setIsSimulating(false);
        }
    };

    const cancelSimulation = async () => {
        if (!simulationId) return;
        
        try {
            await apiClient.cancelSimulation(simulationId);
            setIsSimulating(false);
            addLog('Simulation cancelled');
        } catch (error) {
            console.error('Failed to cancel simulation:', error);
            addLog(`Failed to cancel simulation: ${error.message}`, 'error');
        }
    };
    
    const clearLogs = useCallback(() => {
        setLogs([]);
        addLog('Logs cleared');
    }, [addLog]);
    
    const resetAll = useCallback(() => {
        // Reset to graphene defaults
        setCellVectors([
            [4.260, 0.000, 0.000],
            [0.000, 4.263, 0.000],
            [0.000, 0.000, 6.700]
        ]);
        setBoxRanges([400, 400, 20]);
        setLatticeRanges([
            { min: 0, max: 400 },
            { min: 0, max: 400 },
            { min: 2, max: 3 }
        ]);
        setBasisAtoms([
            { id: generateId(), x: 0.0, y: 0.0, z: 0.0, type: 1 },
            { id: generateId(), x: 0.333, y: 0.667, z: 0.0, type: 2 }
        ]);
        setValidationErrors([]);
        addLog('All parameters reset to defaults');
    }, [addLog]);

    const loadMaterialPreset = useCallback((presetName) => {
        const preset = materialPresets[presetName];
        if (!preset) {
            addLog(`Preset ${presetName} not found`, 'error');
            return;
        }
        
        setCellVectors(preset.cellVectors);
        setBoxRanges(preset.boxRanges);
        setLatticeRanges(preset.latticeRanges);
        setBasisAtoms(preset.basisAtoms.map(atom => ({ ...atom, id: generateId() })));
        addLog(`Loaded ${preset.name} preset`);
    }, [materialPresets, addLog]);

    return (
        <ErrorBoundary>
            <div>
                {/* Header */}
                <header className="header">
                    <div className="header-content">
                        <span className="logo">ARCS</span>
                        <span className="bullet">•</span>
                        <span className="subtitle">Atomic Resolution Collision Simulator Playground</span>
                        <span className="version">version 0.1</span>
                        {healthStatus && (
                            <span className={`health-status ${healthStatus.status}`}>
                                {healthStatus.julia_available ? '🟢' : '🔴'}
                            </span>
                        )}
                    </div>
                    <div className="header-actions">
                        <button className="icon-btn" title="Information" onClick={() => addLog('Application info: DISPLATH v0.1')}>ℹ</button>
                        <button className="icon-btn" title="Validate Parameters" onClick={validateParameters}>✓</button>
                    </div>
                </header>

                {/* Project Bar */}
                <div className="project-bar">
                    <div className="project-info">
                        <div className="project-title-row">
                            <h2 className="project-title">Graphene Sequential Irradiation</h2>
                            <select 
                                className="btn btn-secondary" 
                                onChange={(e) => e.target.value && loadMaterialPreset(e.target.value)}
                                defaultValue=""
                            >
                                <option value="">Load Preset...</option>
                                {Object.entries(materialPresets).map(([key, preset]) => (
                                    <option key={key} value={key}>{preset.name}</option>
                                ))}
                            </select>
                        </div>
                        <div className="project-path">
                            <span>📁</span>
                            <span>/home/Perry/researches/Graphene</span>
                            {simulationId && (
                                <span className="simulation-id">ID: {simulationId.slice(0, 8)}</span>
                            )}
                        </div>
                    </div>
                    <div className="dynamic-load-toggle">
                        <span className="dynamic-load-label">Turn on Dynamic Load</span>
                        <ToggleSwitch 
                            checked={dynamicLoad} 
                            onChange={(e) => setDynamicLoad(e.target.checked)}
                            disabled={isSimulating}
                        />
                    </div>
                </div>

                {/* Validation Errors */}
                {validationErrors.length > 0 && (
                    <div className="validation-errors">
                        <h4>Parameter Validation Errors:</h4>
                        <ul>
                            {validationErrors.map((error, index) => (
                                <li key={index}>{error}</li>
                            ))}
                        </ul>
                    </div>
                )}

                {/* Main Layout */}
                <div className="main-layout">
                    {/* Sidebar */}
                    <aside className={`sidebar ${sidebarCollapsed ? 'collapsed' : ''}`}>
                        <div className="sidebar-toggle">
                            <button className="icon-btn" onClick={() => setSidebarCollapsed(!sidebarCollapsed)}>
                                ☰
                            </button>
                        </div>
                        <div 
                            className={`sidebar-item ${activeTab === 'projectile' ? 'active' : ''}`}
                            onClick={() => setActiveTab('projectile')}
                        >
                            <span className="sidebar-item-icon">🎯</span>
                            <span className="sidebar-item-text">Projectile</span>
                        </div>
                        <div 
                            className={`sidebar-item ${activeTab === 'parameters' ? 'active' : ''}`}
                            onClick={() => setActiveTab('parameters')}
                        >
                            <span className="sidebar-item-icon">⚙️</span>
                            <span className="sidebar-item-text">Parameters</span>
                        </div>
                        <div 
                            className={`sidebar-item ${activeTab === 'process' ? 'active' : ''}`}
                            onClick={() => setActiveTab('process')}
                        >
                            <span className="sidebar-item-icon">🔄</span>
                            <span className="sidebar-item-text">Process</span>
                        </div>
                    </aside>

                    {/* Content Area */}
                    <main className="content-area">
                        {/* Column 1 */}
                        <div className="column">
                            {/* Cell Vectors */}
                            <section className="section">
                                <div className="section-header">
                                    <h3 className="section-title">
                                        Cell Vectors (Å)
                                        <span className="info-icon" title="Unit cell vectors defining crystal structure">ℹ</span>
                                    </h3>
                                </div>
                                <div className="cell-vectors-grid">
                                    {cellVectors.map((row, i) => 
                                        row.map((value, j) => (
                                            <input
                                                key={`cv-${i}-${j}`}
                                                type="number"
                                                className="input-field"
                                                value={value}
                                                onChange={(e) => updateCellVector(i, j, e.target.value)}
                                                step="0.001"
                                                disabled={isSimulating}
                                                placeholder={`a${i+1}[${j+1}]`}
                                            />
                                        ))
                                    )}
                                </div>
                            </section>

                            {/* Basis Atoms */}
                            <section className="section">
                                <div className="section-header">
                                    <h3 className="section-title">Basis Atoms</h3>
                                </div>
                                <div className="basis-atoms-list">
                                    {basisAtoms.map((atom, index) => (
                                        <BasisAtomRow
                                            key={atom.id}
                                            atom={atom}
                                            index={index}
                                            onUpdate={(updated) => updateBasisAtom(atom.id, updated)}
                                            onDelete={() => deleteBasisAtom(atom.id)}
                                        />
                                    ))}
                                </div>
                                <button 
                                    className="btn btn-add" 
                                    onClick={addBasisAtom}
                                    disabled={isSimulating}
                                >
                                    + New Basis Atom
                                </button>
                            </section>

                            {/* Element Library */}
                            <section className="section">
                                <div className="section-header">
                                    <h3 className="section-title">Element Library</h3>
                                </div>
                                <div className="element-library">
                                    {elements.map(element => (
                                        <ElementCard
                                            key={element.number}
                                            number={element.number}
                                            symbol={element.symbol}
                                            properties={element.properties}
                                        />
                                    ))}
                                </div>
                                <button className="btn btn-add" disabled={isSimulating}>
                                    + New Element
                                </button>
                            </section>
                        </div>

                        {/* Column 2 */}
                        <div className="column">
                            {/* Box Ranges */}
                            <section className="section">
                                <div className="section-header">
                                    <h3 className="section-title">
                                        Box Ranges
                                        <span className="info-icon" title="Number of unit cells in each direction">ℹ</span>
                                    </h3>
                                </div>
                                <div className="box-ranges-inputs">
                                    {boxRanges.map((value, index) => (
                                        <input
                                            key={`br-${index}`}
                                            type="number"
                                            className="input-field"
                                            value={value}
                                            onChange={(e) => updateBoxRange(index, e.target.value)}
                                            disabled={isSimulating}
                                            min="1"
                                            placeholder={['X', 'Y', 'Z'][index]}
                                        />
                                    ))}
                                </div>
                                <div className="dimensions-info">
                                    <div>Dimensions: {dimensions.x} Å × {dimensions.y} Å × {dimensions.z} Å</div>
                                    <div>Center Point: ({centerPoint.x}, {centerPoint.y}, {centerPoint.z})</div>
                                </div>
                            </section>

                            {/* Lattice Range */}
                            <section className="section">
                                <div className="section-header">
                                    <h3 className="section-title">
                                        Lattice Range
                                        <span className="info-icon" title="Range of lattice cells to include">ℹ</span>
                                    </h3>
                                </div>
                                {latticeRanges.map((range, index) => (
                                    <div key={`lr-${index}`} className="lattice-range-row">
                                        <input
                                            type="number"
                                            className="input-field"
                                            value={range.min}
                                            onChange={(e) => updateLatticeRange(index, 'min', e.target.value)}
                                            disabled={isSimulating}
                                            placeholder="Min"
                                        />
                                        <span className="lattice-range-label">to</span>
                                        <input
                                            type="number"
                                            className="input-field"
                                            value={range.max}
                                            onChange={(e) => updateLatticeRange(index, 'max', e.target.value)}
                                            disabled={isSimulating}
                                            placeholder="Max"
                                        />
                                        <span className="lattice-range-value">
                                            {(range.min * (cellVectors[index]?.[index] || 1)).toFixed(1)} Å
                                        </span>
                                        <span className="lattice-range-value">
                                            {(range.max * (cellVectors[index]?.[index] || 1)).toFixed(1)} Å
                                        </span>
                                    </div>
                                ))}
                            </section>

                            {/* Simulation Progress */}
                            {simulationStatus && (
                                <section className="section">
                                    <div className="section-header">
                                        <h3 className="section-title">Simulation Progress</h3>
                                    </div>
                                    <div className="simulation-progress">
                                        <div className="progress-bar">
                                            <div 
                                                className="progress-fill" 
                                                style={{ width: `${simulationStatus.progress || 0}%` }}
                                            ></div>
                                        </div>
                                        <div className="progress-info">
                                            <span>Status: {simulationStatus.status}</span>
                                            <span>Progress: {simulationStatus.progress || 0}%</span>
                                            {simulationStatus.results && (
                                                <span>Vacancies: {simulationStatus.results.totalVacancies}</span>
                                            )}
                                        </div>
                                    </div>
                                </section>
                            )}
                        </div>

                        {/* Column 3 */}
                        <div className="column">
                            {/* Crystal Visualization */}
                            <section className="section">
                                <div className="section-header">
                                    <h3 className="section-title">
                                        Crystal Visualization
                                        <span className="info-icon" title="3D view of crystal structure">ℹ</span>
                                    </h3>
                                    <button className="btn btn-secondary" onClick={resetAll} disabled={isSimulating}>
                                        Reset
                                    </button>
                                </div>
                                <div className="crystal-viz">
                                    <CrystalVisualization 
                                        cellVectors={cellVectors}
                                        basisAtoms={basisAtoms}
                                    />
                                </div>
                            </section>

                            {/* System Logs */}
                            <section className="section">
                                <div className="section-header">
                                    <h3 className="section-title">System Logs</h3>
                                </div>
                                <div className="system-logs">
                                    {logs.slice(-50).map((log, index) => (
                                        <div key={index} className={`log-entry log-${log.level || 'info'}`}>
                                            <span className="log-time">[{log.time}]</span> {log.message}
                                        </div>
                                    ))}
                                </div>
                            </section>

                            {/* Control Buttons */}
                            <div className="control-buttons">
                                <button className="btn btn-secondary" onClick={clearLogs}>
                                    Clear
                                </button>
                                {isSimulating ? (
                                    <button className="btn btn-secondary" onClick={cancelSimulation}>
                                        ⏹ Cancel
                                    </button>
                                ) : (
                                    <button 
                                        className="btn btn-primary" 
                                        onClick={runSimulation}
                                        disabled={!healthStatus?.julia_available}
                                    >
                                        ▶ Run Simulation
                                    </button>
                                )}
                                <button className="btn btn-secondary" onClick={resetAll} disabled={isSimulating}>
                                    ↻ Reset
                                </button>
                            </div>
                        </div>
                    </main>
                </div>
            </div>
        </ErrorBoundary>
    );
};

// Create React root and render the app
const container = document.getElementById('root');
const root = ReactDOM.createRoot ? ReactDOM.createRoot(container) : null;

if (root) {
    root.render(<App />);
} else {
    // Fallback for older React versions
    ReactDOM.render(<App />, container);
}