console.log('DISPLATH GUI loaded');

// DISPLATH BCA Simulator - Modern JavaScript

// Global State
const state = {
    currentTab: 'material',
    simulation: {
        isRunning: false,
        progress: 0,
        stats: {
            totalCollisions: 0,
            vacancies: 0,
            displacements: 0,
            interstitials: 0,
            sputteredAtoms: 0,
            cascadeDepth: 0
        }
    },
    parameters: {
        material: {
            basisAtoms: [
                { position: [0.000, 0.000, 0.000], element: 1 },
                { position: [0.500, 0.500, 0.000], element: 1 },
                { position: [0.500, 0.000, 0.500], element: 1 },
                { position: [0.000, 0.500, 0.500], element: 1 },
                { position: [0.250, 0.250, 0.250], element: 1 },
                { position: [0.750, 0.750, 0.250], element: 1 },
                { position: [0.750, 0.250, 0.750], element: 1 },
                { position: [0.250, 0.750, 0.750], element: 1 }
            ],
            elements: {
                1: { name: 'Si', radius: 2.20, mass: 28.09, Z: 14, dte: 35.0, bde: 35.0 },
                2: { name: 'B', radius: 1.00, mass: 10.81, Z: 5, dte: 25.0, bde: 25.0 },
                3: { name: 'Ar', radius: 1.88, mass: 39.95, Z: 18, dte: 0.1, bde: 0.1 }
            }
        }
    }
};

// Initialize on DOM load
document.addEventListener('DOMContentLoaded', function() {
    initializeTabs();
    initializeMatrixInputs();
    initializeBasisAtoms();
    initializeElementLibrary();
    initializeDirectionCalculator();
    initializeAdvancedSection();
    initializeActionButtons();
    updateAtomCount();
    addLog('System initialized. Ready for simulation.');
});

// Tab Management
function initializeTabs() {
    const tabBtns = document.querySelectorAll('.tab-btn');
    const tabPanes = document.querySelectorAll('.tab-pane');

    tabBtns.forEach(btn => {
        btn.addEventListener('click', function() {
            const targetTab = this.dataset.tab;
            
            tabBtns.forEach(b => b.classList.remove('active'));
            this.classList.add('active');
            
            tabPanes.forEach(pane => pane.classList.remove('active'));
            document.getElementById(targetTab).classList.add('active');
            
            state.currentTab = targetTab;
        });
    });
}

// Matrix Input Handling
function initializeMatrixInputs() {
    const primaryInputs = ['a1x', 'a1y', 'a1z', 'a2x', 'a2y', 'a2z', 'a3x', 'a3y', 'a3z'];
    primaryInputs.forEach(id => {
        const input = document.getElementById(id);
        if (input) {
            input.addEventListener('input', updateMatrixInfo);
        }
    });

    const rangeInputs = ['xMin', 'xMax', 'yMin', 'yMax', 'zMin', 'zMax'];
    rangeInputs.forEach(id => {
        const input = document.getElementById(id);
        if (input) {
            input.addEventListener('input', updateAtomCount);
        }
    });

    updateMatrixInfo();
}

function updateMatrixInfo() {
    const matrix = [
        [parseFloat(document.getElementById('a1x').value) || 0, parseFloat(document.getElementById('a1y').value) || 0, parseFloat(document.getElementById('a1z').value) || 0],
        [parseFloat(document.getElementById('a2x').value) || 0, parseFloat(document.getElementById('a2y').value) || 0, parseFloat(document.getElementById('a2z').value) || 0],
        [parseFloat(document.getElementById('a3x').value) || 0, parseFloat(document.getElementById('a3y').value) || 0, parseFloat(document.getElementById('a3z').value) || 0]
    ];

    const det = calculateDeterminant(matrix);
    const volume = Math.abs(det);
    const isOrthogonal = isMatrixOrthogonal(matrix);
    
    document.getElementById('matrixDet').textContent = det.toFixed(3);
    document.getElementById('unitCellVol').textContent = volume.toFixed(1);
    
    const orthogonalIndicator = document.getElementById('orthogonalIndicator');
    if (isOrthogonal) {
        orthogonalIndicator.textContent = '✓ Orthogonal';
        orthogonalIndicator.style.color = 'var(--color-success)';
    } else {
        orthogonalIndicator.textContent = '✗ Non-orthogonal';
        orthogonalIndicator.style.color = 'var(--color-warning)';
    }

    updateAtomCount();
}

function calculateDeterminant(matrix) {
    const [[a, b, c], [d, e, f], [g, h, i]] = matrix;
    return a * (e * i - f * h) - b * (d * i - f * g) + c * (d * h - e * g);
}

function isMatrixOrthogonal(matrix) {
    const [[a1, a2, a3], [b1, b2, b3], [c1, c2, c3]] = matrix;
    const threshold = 1e-6;
    return Math.abs(a2) < threshold && Math.abs(a3) < threshold &&
           Math.abs(b1) < threshold && Math.abs(b3) < threshold &&
           Math.abs(c1) < threshold && Math.abs(c2) < threshold;
}

function updateAtomCount() {
    const xRange = (parseInt(document.getElementById('xMax').value) || 0) - (parseInt(document.getElementById('xMin').value) || 0);
    const yRange = (parseInt(document.getElementById('yMax').value) || 0) - (parseInt(document.getElementById('yMin').value) || 0);
    const zRange = (parseInt(document.getElementById('zMax').value) || 0) - (parseInt(document.getElementById('zMin').value) || 0);
    
    const basisCount = state.parameters.material.basisAtoms.length;
    const totalAtoms = xRange * yRange * zRange * basisCount;
    
    let displayText;
    if (totalAtoms > 1e9) {
        displayText = `~${(totalAtoms / 1e9).toFixed(1)}B`;
    } else if (totalAtoms > 1e6) {
        displayText = `~${(totalAtoms / 1e6).toFixed(1)}M`;
    } else if (totalAtoms > 1e3) {
        displayText = `~${(totalAtoms / 1e3).toFixed(1)}K`;
    } else {
        displayText = totalAtoms.toString();
    }
    
    document.getElementById('atomCount').textContent = displayText;
}

// Basis Atoms Management
function initializeBasisAtoms() {
    renderBasisAtoms();
    document.getElementById('addBasisAtom').addEventListener('click', addBasisAtom);
}

function renderBasisAtoms() {
    const container = document.getElementById('basisAtomsList');
    container.innerHTML = '';

    state.parameters.material.basisAtoms.forEach((atom, index) => {
        const atomDiv = document.createElement('div');
        atomDiv.className = 'basis-atom-item';
        atomDiv.innerHTML = `
            <div class="atom-info">
                <span class="atom-label">Atom ${index + 1}:</span>
                <div class="atom-position">
                    <input type="number" step="0.001" value="${atom.position[0]}" onchange="updateBasisAtom(${index}, 0, this.value)">
                    <input type="number" step="0.001" value="${atom.position[1]}" onchange="updateBasisAtom(${index}, 1, this.value)">
                    <input type="number" step="0.001" value="${atom.position[2]}" onchange="updateBasisAtom(${index}, 2, this.value)">
                </div>
                <select onchange="updateBasisAtomElement(${index}, this.value)">
                    ${Object.entries(state.parameters.material.elements).map(([id, el]) => 
                        `<option value="${id}" ${atom.element == id ? 'selected' : ''}>${el.name}</option>`
                    ).join('')}
                </select>
                <button class="btn-icon" onclick="removeBasisAtom(${index})" title="Remove">✕</button>
            </div>
        `;
        container.appendChild(atomDiv);
    });
}

function addBasisAtom() {
    state.parameters.material.basisAtoms.push({
        position: [0.0, 0.0, 0.0],
        element: 1
    });
    renderBasisAtoms();
    updateAtomCount();
}

function updateBasisAtom(index, coord, value) {
    state.parameters.material.basisAtoms[index].position[coord] = parseFloat(value);
    updateAtomCount();
}

function updateBasisAtomElement(index, elementId) {
    state.parameters.material.basisAtoms[index].element = parseInt(elementId);
}

function removeBasisAtom(index) {
    state.parameters.material.basisAtoms.splice(index, 1);
    renderBasisAtoms();
    updateAtomCount();
}

// Element Library
function initializeElementLibrary() {
    renderElementLibrary();
    document.getElementById('addElement').addEventListener('click', addElement);
}

function renderElementLibrary() {
    const container = document.getElementById('elementLibrary');
    container.innerHTML = '';

    Object.entries(state.parameters.material.elements).forEach(([id, element]) => {
        const elementDiv = document.createElement('div');
        elementDiv.className = 'element-item';
        elementDiv.innerHTML = `
            <div class="element-info">
                <span class="element-name">${element.name}</span>
                <div class="element-properties">
                    <span>R: ${element.radius}Å</span>
                    <span>M: ${element.mass}u</span>
                    <span>Z: ${element.Z}</span>
                    <span>DTE: ${element.dte}eV</span>
                </div>
                <button class="btn-icon" onclick="editElement(${id})" title="Edit">⚙️</button>
            </div>
        `;
        container.appendChild(elementDiv);
    });
}

function addElement() {
    const newId = Math.max(...Object.keys(state.parameters.material.elements).map(Number)) + 1;
    state.parameters.material.elements[newId] = {
        name: 'New',
        radius: 1.0,
        mass: 1.0,
        Z: 1,
        dte: 1.0,
        bde: 1.0
    };
    renderElementLibrary();
    renderBasisAtoms();
}

function editElement(id) {
    alert(`Edit element ${id} - Would open element editor modal`);
}

// Direction Calculator
function initializeDirectionCalculator() {
    const thetaInput = document.getElementById('dirTheta');
    const phiInput = document.getElementById('dirPhi');
    
    if (thetaInput && phiInput) {
        thetaInput.addEventListener('input', updateDirectionVector);
        phiInput.addEventListener('input', updateDirectionVector);
        updateDirectionVector();
    }
}

function updateDirectionVector() {
    const theta = parseFloat(document.getElementById('dirTheta').value) || 0;
    const phi = parseFloat(document.getElementById('dirPhi').value) || 0;
    
    const thetaRad = theta * Math.PI / 180;
    const phiRad = phi * Math.PI / 180;
    
    const x = Math.sin(thetaRad) * Math.cos(phiRad);
    const y = Math.sin(thetaRad) * Math.sin(phiRad);
    const z = -Math.cos(thetaRad);
    
    document.getElementById('dirVecX').textContent = x.toFixed(3);
    document.getElementById('dirVecY').textContent = y.toFixed(3);
    document.getElementById('dirVecZ').textContent = z.toFixed(3);
}

// Advanced Section
function initializeAdvancedSection() {
    const toggle = document.getElementById('advancedToggle');
    const content = document.getElementById('advancedContent');
    
    toggle.addEventListener('click', function() {
        const isExpanded = this.classList.contains('expanded');
        
        if (isExpanded) {
            this.classList.remove('expanded');
            content.classList.remove('expanded');
        } else {
            this.classList.add('expanded');
            content.classList.add('expanded');
        }
    });
}

// Action Buttons
function initializeActionButtons() {
    document.getElementById('validateParams').addEventListener('click', validateParameters);
    document.getElementById('runSimulation').addEventListener('click', runSimulation);
    document.getElementById('resetParams').addEventListener('click', resetParameters);
    document.getElementById('randomizeSeed').addEventListener('click', randomizeSeed);
    document.getElementById('clearLogs').addEventListener('click', clearLogs);
    document.getElementById('exportLogs').addEventListener('click', exportLogs);
}

function validateParameters() {
    addLog('Validating parameters...');
    const errors = [];
    
    const det = parseFloat(document.getElementById('matrixDet').textContent);
    if (Math.abs(det) < 1e-10) {
        errors.push('Primary vectors matrix is singular');
    }
    
    if (errors.length === 0) {
        addLog('✓ Parameters validation passed', 'success');
    } else {
        errors.forEach(error => addLog(`✗ ${error}`, 'error'));
    }
    
    return errors.length === 0;
}

async function runSimulation() {
    if (!validateParameters()) {
        addLog('Cannot run simulation: validation failed', 'error');
        return;
    }
    
    if (state.simulation.isRunning) {
        addLog('Simulation is already running', 'warning');
        return;
    }
    
    try {
        state.simulation.isRunning = true;
        updateSimulationButtons();
        
        addLog('Starting simulation...');
        
        // Simulate progress
        for (let i = 0; i <= 100; i += 10) {
            await new Promise(resolve => setTimeout(resolve, 500));
            updateProgress(i);
            updateStatistics({
                totalCollisions: Math.floor(Math.random() * 50000),
                vacancies: Math.floor(Math.random() * 500)
            });
        }
        
        addLog('✓ Simulation completed successfully', 'success');
        
    } catch (error) {
        addLog(`✗ Simulation failed: ${error.message}`, 'error');
    } finally {
        state.simulation.isRunning = false;
        updateSimulationButtons();
    }
}

function resetParameters() {
    addLog('Resetting parameters...');
    document.getElementById('a1x').value = 5.431;
    document.getElementById('a1y').value = 0.000;
    document.getElementById('a1z').value = 0.000;
    updateMatrixInfo();
    addLog('Parameters reset');
}

function randomizeSeed() {
    const newSeed = Math.floor(Math.random() * 10000);
    document.getElementById('randomSeed').value = newSeed;
    addLog(`Random seed set to: ${newSeed}`);
}

// Progress and Statistics
function updateProgress(percent) {
    state.simulation.progress = percent;
    document.getElementById('progressFill').style.width = `${percent}%`;
    document.getElementById('progressPercent').textContent = percent.toFixed(1);
}

function updateStatistics(stats) {
    if (stats) {
        Object.assign(state.simulation.stats, stats);
    }
    
    const s = state.simulation.stats;
    document.getElementById('totalCollisions').textContent = s.totalCollisions.toLocaleString();
    document.getElementById('vacancies').textContent = s.vacancies.toLocaleString();
    document.getElementById('displacements').textContent = s.displacements.toLocaleString();
    document.getElementById('interstitials').textContent = s.interstitials.toLocaleString();
    document.getElementById('sputteredAtoms').textContent = s.sputteredAtoms.toLocaleString();
    document.getElementById('cascadeDepth').textContent = `${s.cascadeDepth.toFixed(1)} Å`;
}

function updateSimulationButtons() {
    const runBtn = document.getElementById('runSimulation');
    
    if (state.simulation.isRunning) {
        runBtn.innerHTML = '<span>⏸</span> Running...';
        runBtn.disabled = true;
    } else {
        runBtn.innerHTML = '<span>▶</span> Run Simulation';
        runBtn.disabled = false;
    }
}

// Logging
function addLog(message, type = 'info') {
    const logsContent = document.getElementById('logsContent');
    const time = new Date().toLocaleTimeString();
    
    const logEntry = document.createElement('div');
    logEntry.className = 'log-entry';
    logEntry.innerHTML = `
        <span class="log-time">[${time}]</span>
        <span class="log-message ${type}">${message}</span>
    `;
    
    logsContent.appendChild(logEntry);
    logsContent.scrollTop = logsContent.scrollHeight;
}

function clearLogs() {
    document.getElementById('logsContent').innerHTML = '';
    addLog('Logs cleared');
}

function exportLogs() {
    const logs = document.getElementById('logsContent').innerText;
    const blob = new Blob([logs], { type: 'text/plain' });
    const url = URL.createObjectURL(blob);
    
    const a = document.createElement('a');
    a.href = url;
    a.download = `displath_logs_${new Date().toISOString().slice(0, 19).replace(/:/g, '-')}.txt`;
    a.click();
    
    URL.revokeObjectURL(url);
    addLog('Logs exported');
}

// Global functions for inline handlers
window.updateBasisAtom = updateBasisAtom;
window.updateBasisAtomElement = updateBasisAtomElement;
window.removeBasisAtom = removeBasisAtom;
window.editElement = editElement;
