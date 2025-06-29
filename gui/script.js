// Global state
let currentChart = null;

// DOM Elements
const elements = {
    // Tabs
    tabButtons: document.querySelectorAll('.tab-button'),
    tabPanes: document.querySelectorAll('.tab-pane'),
    
    // Controls
    presetSelect: document.getElementById('presetSelect'),
    runSimulation: document.getElementById('runSimulation'),
    resetForm: document.getElementById('resetForm'),
    isDynamicLoad: document.getElementById('isDynamicLoad'),
    dynamicLoadOptions: document.getElementById('dynamicLoadOptions'),
    
    // Basis and Element management
    addBasisAtom: document.getElementById('addBasisAtom'),
    addElementType: document.getElementById('addElementType'),
    basisContainer: document.getElementById('basisContainer'),
    typeContainer: document.getElementById('typeContainer'),
    
    // Results
    loadingIndicator: document.getElementById('loadingIndicator'),
    errorDisplay: document.getElementById('errorDisplay'),
    resultsDisplay: document.getElementById('resultsDisplay'),
    errorMessage: document.getElementById('errorMessage'),
    
    // Stats
    totalRuns: document.getElementById('totalRuns'),
    meanVacancies: document.getElementById('meanVacancies'),
    stdVacancies: document.getElementById('stdVacancies'),
    rawDataTable: document.getElementById('rawDataTable'),
    histogramChart: document.getElementById('histogramChart')
};

// Initialize the application
document.addEventListener('DOMContentLoaded', function() {
    initializeTabs();
    initializeEventListeners();
    initializeFormHandlers();
    loadPresets();
});

// Tab Management
function initializeTabs() {
    elements.tabButtons.forEach(button => {
        button.addEventListener('click', function() {
            const tabId = this.dataset.tab;
            switchTab(tabId);
        });
    });
}

function switchTab(tabId) {
    // Update tab buttons
    elements.tabButtons.forEach(btn => btn.classList.remove('active'));
    document.querySelector(`[data-tab="${tabId}"]`).classList.add('active');
    
    // Update tab panes
    elements.tabPanes.forEach(pane => pane.classList.remove('active'));
    document.getElementById(tabId).classList.add('active');
}

// Event Listeners
function initializeEventListeners() {
    elements.presetSelect.addEventListener('change', handlePresetChange);
    elements.runSimulation.addEventListener('click', handleRunSimulation);
    elements.resetForm.addEventListener('click', handleResetForm);
    elements.isDynamicLoad.addEventListener('change', handleDynamicLoadToggle);
    elements.addBasisAtom.addEventListener('click', addBasisAtom);
    elements.addElementType.addEventListener('click', addElementType);
}

function initializeFormHandlers() {
    // Auto-update ion type options when element types change
    document.addEventListener('change', function(e) {
        if (e.target.classList.contains('element-name')) {
            updateIonTypeOptions();
        }
    });
}

// Preset Management
async function loadPresets() {
    try {
        const response = await fetch('/api/presets');
        const presets = await response.json();
        
        Object.keys(presets).forEach(key => {
            const option = document.createElement('option');
            option.value = key;
            option.textContent = presets[key].name;
            elements.presetSelect.appendChild(option);
        });
    } catch (error) {
        console.error('Failed to load presets:', error);
    }
}

function handlePresetChange() {
    const presetKey = elements.presetSelect.value;
    if (!presetKey) return;
    
    fetch('/api/presets')
        .then(response => response.json())
        .then(presets => {
            const preset = presets[presetKey];
            if (preset) {
                applyPreset(preset);
            }
        })
        .catch(error => {
            console.error('Failed to apply preset:', error);
            showError('Failed to load preset configuration');
        });
}

function applyPreset(preset) {
    // Primary vectors
    const primaryVectorInputs = document.querySelectorAll('#primaryVectors input');
    preset.primaryVectors.forEach((value, index) => {
        if (primaryVectorInputs[index]) {
            primaryVectorInputs[index].value = value;
        }
    });
    
    // Lattice ranges
    const ranges = ['xMin', 'xMax', 'yMin', 'yMax', 'zMin', 'zMax'];
    preset.latticeRanges.forEach((value, index) => {
        const element = document.getElementById(ranges[index]);
        if (element) element.value = value;
    });
    
    // Clear existing basis atoms and element types
    elements.basisContainer.innerHTML = '';
    elements.typeContainer.innerHTML = '';
    
    // Add basis atoms
    const basisCount = preset.basis.length / 3;
    for (let i = 0; i < basisCount; i++) {
        const x = preset.basis[i * 3];
        const y = preset.basis[i * 3 + 1];
        const z = preset.basis[i * 3 + 2];
        const type = preset.basisTypes[i];
        addBasisAtom(x, y, z, type);
    }
    
    // Add element types
    Object.keys(preset.typeDict).forEach(typeId => {
        const element = preset.typeDict[typeId];
        addElementType(element.name, element.dte, element.bde, typeId);
    });
    
    updateIonTypeOptions();
}

// Dynamic Load Toggle
function handleDynamicLoadToggle() {
    if (elements.isDynamicLoad.checked) {
        elements.dynamicLoadOptions.style.display = 'block';
    } else {
        elements.dynamicLoadOptions.style.display = 'none';
    }
}

// Basis Atom Management
function addBasisAtom(x = 0, y = 0, z = 0, type = 1) {
    const atomDiv = document.createElement('div');
    atomDiv.className = 'basis-atom';
    
    const atomIndex = elements.basisContainer.children.length + 1;
    
    atomDiv.innerHTML = `
        <span>Atom ${atomIndex}:</span>
        <input type="number" step="0.001" value="${x}" class="basis-x">
        <input type="number" step="0.001" value="${y}" class="basis-y">
        <input type="number" step="0.001" value="${z}" class="basis-z">
        <select class="basis-type">
            <option value="1" ${type == 1 ? 'selected' : ''}>Type 1</option>
            <option value="2" ${type == 2 ? 'selected' : ''}>Type 2</option>
        </select>
        <button type="button" class="btn-remove" onclick="removeBasisAtom(this)">
            <i class="fas fa-trash"></i>
        </button>
    `;
    
    elements.basisContainer.appendChild(atomDiv);
    updateBasisAtomLabels();
}

function removeBasisAtom(button) {
    button.parentElement.remove();
    updateBasisAtomLabels();
}

function updateBasisAtomLabels() {
    const atoms = elements.basisContainer.querySelectorAll('.basis-atom');
    atoms.forEach((atom, index) => {
        atom.querySelector('span').textContent = `Atom ${index + 1}:`;
    });
}

// Element Type Management
function addElementType(name = 'C', dte = 19.96, bde = 19.96, typeId = null) {
    const typeDiv = document.createElement('div');
    typeDiv.className = 'element-type';
    
    const typeIndex = typeId || (elements.typeContainer.children.length + 1);
    
    typeDiv.innerHTML = `
        <span>Type ${typeIndex}:</span>
        <input type="text" value="${name}" class="element-name" placeholder="Symbol">
        <input type="number" step="0.01" value="${dte}" class="element-dte" placeholder="DTE (eV)">
        <input type="number" step="0.01" value="${bde}" class="element-bde" placeholder="BDE (eV)">
        <button type="button" class="btn-remove" onclick="removeElementType(this)">
            <i class="fas fa-trash"></i>
        </button>
    `;
    
    typeDiv.dataset.typeId = typeIndex;
    elements.typeContainer.appendChild(typeDiv);
    updateElementTypeLabels();
    updateIonTypeOptions();
}

function removeElementType(button) {
    button.parentElement.remove();
    updateElementTypeLabels();
    updateIonTypeOptions();
}

function updateElementTypeLabels() {
    const types = elements.typeContainer.querySelectorAll('.element-type');
    types.forEach((type, index) => {
        const typeId = index + 1;
        type.querySelector('span').textContent = `Type ${typeId}:`;
        type.dataset.typeId = typeId;
    });
    
    // Update basis type options
    const basisSelects = document.querySelectorAll('.basis-type');
    basisSelects.forEach(select => {
        const currentValue = select.value;
        select.innerHTML = '';
        
        types.forEach((type, index) => {
            const option = document.createElement('option');
            option.value = index + 1;
            option.textContent = `Type ${index + 1}`;
            if (option.value === currentValue) option.selected = true;
            select.appendChild(option);
        });
    });
}

function updateIonTypeOptions() {
    const ionTypeSelect = document.getElementById('ionType');
    const currentValue = ionTypeSelect.value;
    ionTypeSelect.innerHTML = '';
    
    const types = elements.typeContainer.querySelectorAll('.element-type');
    types.forEach((type, index) => {
        const typeId = index + 1;
        const name = type.querySelector('.element-name').value;
        const option = document.createElement('option');
        option.value = typeId;
        option.textContent = `Type ${typeId} (${name})`;
        if (option.value === currentValue) option.selected = true;
        ionTypeSelect.appendChild(option);
    });
}

// Form Handling
function handleResetForm() {
    if (confirm('Are you sure you want to reset all form data?')) {
        location.reload();
    }
}

function collectFormData() {
    // Primary vectors
    const primaryVectorInputs = document.querySelectorAll('#primaryVectors input');
    const primaryVectors = Array.from(primaryVectorInputs).map(input => parseFloat(input.value));
    
    // Lattice ranges
    const latticeRanges = [
        parseInt(document.getElementById('xMin').value),
        parseInt(document.getElementById('xMax').value),
        parseInt(document.getElementById('yMin').value),
        parseInt(document.getElementById('yMax').value),
        parseInt(document.getElementById('zMin').value),
        parseInt(document.getElementById('zMax').value)
    ];
    
    // Box sizes
    const boxSizes = [
        parseInt(document.getElementById('boxX').value),
        parseInt(document.getElementById('boxY').value),
        parseInt(document.getElementById('boxZ').value)
    ];
    
    // Input grid vectors (same as primary vectors for now)
    const inputGridVectors = [...primaryVectors];
    
    // Basis atoms
    const basisAtoms = Array.from(elements.basisContainer.querySelectorAll('.basis-atom'));
    const basis = [];
    const basisTypes = [];
    
    basisAtoms.forEach(atom => {
        basis.push(
            parseFloat(atom.querySelector('.basis-x').value),
            parseFloat(atom.querySelector('.basis-y').value),
            parseFloat(atom.querySelector('.basis-z').value)
        );
        basisTypes.push(parseInt(atom.querySelector('.basis-type').value));
    });
    
    // Element types
    const elementTypes = Array.from(elements.typeContainer.querySelectorAll('.element-type'));
    const typeDict = {};
    
    elementTypes.forEach((type, index) => {
        const typeId = index + 1;
        typeDict[typeId] = {
            name: type.querySelector('.element-name').value,
            dte: parseFloat(type.querySelector('.element-dte').value),
            bde: parseFloat(type.querySelector('.element-bde').value)
        };
    });
    
    // Ion parameters
    const ionDirection = [
        parseFloat(document.getElementById('dirX').value),
        parseFloat(document.getElementById('dirY').value),
        parseFloat(document.getElementById('dirZ').value)
    ];
    
    const ionPosition = [
        parseFloat(document.getElementById('posX').value),
        parseFloat(document.getElementById('posY').value),
        parseFloat(document.getElementById('posZ').value)
    ];
    
    return {
        primaryVectors,
        latticeRanges,
        basisTypes,
        basis,
        boxSizes,
        inputGridVectors,
        pMax: parseFloat(document.getElementById('pMax').value),
        vacancyRecoverDistance: parseFloat(document.getElementById('vacancyRecoverDistance').value),
        temperature: parseFloat(document.getElementById('temperature').value),
        DebyeTemperature: parseFloat(document.getElementById('DebyeTemperature').value),
        stopEnergy: parseFloat(document.getElementById('stopEnergy').value),
        ionType: parseInt(document.getElementById('ionType').value),
        ionEnergy: parseFloat(document.getElementById('ionEnergy').value),
        ionDirection,
        ionPosition,
        nRuns: parseInt(document.getElementById('nRuns').value),
        isDynamicLoad: elements.isDynamicLoad.checked,
        nCascadeEveryLoad: parseInt(document.getElementById('nCascadeEveryLoad').value),
        outputName: document.getElementById('outputName').value,
        isDumpInCascade: document.getElementById('isDumpInCascade').checked,
        typeDict
    };
}

// Simulation Handling
async function handleRunSimulation() {
    try {
        showLoading();
        
        const formData = collectFormData();
        
        // Validate form data
        if (!validateFormData(formData)) {
            throw new Error('Please check your input parameters');
        }
        
        const response = await fetch('/api/simulate', {
            method: 'POST',
            headers: {
                'Content-Type': 'application/json'
            },
            body: JSON.stringify(formData)
        });
        
        const result = await response.json();
        
        if (result.success) {
            showResults(result.results);
        } else {
            throw new Error(result.error || 'Simulation failed');
        }
        
    } catch (error) {
        console.error('Simulation error:', error);
        showError(error.message);
    } finally {
        hideLoading();
    }
}

function validateFormData(data) {
    // Basic validation
    if (data.nRuns <= 0 || data.nRuns > 10000) {
        alert('Number of runs must be between 1 and 10000');
        return false;
    }
    
    if (data.ionEnergy <= 0) {
        alert('Ion energy must be positive');
        return false;
    }
    
    if (data.basis.length === 0) {
        alert('At least one basis atom is required');
        return false;
    }
    
    if (Object.keys(data.typeDict).length === 0) {
        alert('At least one element type is required');
        return false;
    }
    
    return true;
}

// Results Display
function showResults(results) {
    hideError();
    
    // Update stats
    elements.totalRuns.textContent = results.totalRuns;
    elements.meanVacancies.textContent = results.meanVacancies.toFixed(2);
    elements.stdVacancies.textContent = results.stdVacancies.toFixed(2);
    
    // Create histogram
    createHistogram(results.allVacancies);
    
    // Create data table
    createDataTable(results.allVacancies);
    
    elements.resultsDisplay.style.display = 'block';
}

function createHistogram(vacancies) {
    const ctx = elements.histogramChart.getContext('2d');
    
    // Destroy existing chart
    if (currentChart) {
        currentChart.destroy();
    }
    
    // Calculate histogram data
    const maxVacancies = Math.max(...vacancies);
    const minVacancies = Math.min(...vacancies);
    const binCount = Math.min(20, maxVacancies - minVacancies + 1);
    
    const bins = Array(binCount).fill(0);
    const binWidth = (maxVacancies - minVacancies) / (binCount - 1);
    
    vacancies.forEach(v => {
        const binIndex = Math.min(Math.floor((v - minVacancies) / binWidth), binCount - 1);
        bins[binIndex]++;
    });
    
    const labels = bins.map((_, i) => {
        const binStart = minVacancies + i * binWidth;
        return Math.round(binStart);
    });
    
    currentChart = new Chart(ctx, {
        type: 'bar',
        data: {
            labels: labels,
            datasets: [{
                label: 'Frequency',
                data: bins,
                backgroundColor: 'rgba(37, 99, 235, 0.6)',
                borderColor: 'rgba(37, 99, 235, 1)',
                borderWidth: 1
            }]
        },
        options: {
            responsive: true,
            plugins: {
                title: {
                    display: true,
                    text: 'Vacancy Distribution'
                }
            },
            scales: {
                x: {
                    title: {
                        display: true,
                        text: 'Number of Vacancies'
                    }
                },
                y: {
                    title: {
                        display: true,
                        text: 'Frequency'
                    },
                    beginAtZero: true
                }
            }
        }
    });
}

function createDataTable(vacancies) {
    const table = document.createElement('table');
    table.innerHTML = `
        <thead>
            <tr>
                <th>Run #</th>
                <th>Vacancies</th>
            </tr>
        </thead>
        <tbody>
            ${vacancies.map((v, i) => `
                <tr>
                    <td>${i + 1}</td>
                    <td>${v}</td>
                </tr>
            `).join('')}
        </tbody>
    `;
    
    elements.rawDataTable.innerHTML = '';
    elements.rawDataTable.appendChild(table);
}

// UI State Management
function showLoading() {
    elements.loadingIndicator.style.display = 'block';
    elements.errorDisplay.style.display = 'none';
    elements.resultsDisplay.style.display = 'none';
    elements.runSimulation.disabled = true;
}

function hideLoading() {
    elements.loadingIndicator.style.display = 'none';
    elements.runSimulation.disabled = false;
}

function showError(message) {
    elements.errorMessage.textContent = message;
    elements.errorDisplay.style.display = 'flex';
    elements.resultsDisplay.style.display = 'none';
}

function hideError() {
    elements.errorDisplay.style.display = 'none';
}

// Initialize with default data
setTimeout(() => {
    // Add default basis atoms for graphene
    addBasisAtom(0, 0, 0, 1);
    addBasisAtom(1/3, 0, 0, 1);
    addBasisAtom(1/2, 1/2, 0, 1);
    addBasisAtom(5/6, 1/2, 0, 1);
    
    // Add default element types
    addElementType('C', 19.96, 19.96, 1);
    addElementType('Ne', 1.0, 1.0, 2);
}, 100);