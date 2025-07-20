#!/usr/bin/env python3
import os
import sys
import subprocess
import threading
import queue
import json
from flask import Flask, render_template, request, jsonify
from flask_socketio import SocketIO, emit, join_room
import time
from pathlib import Path
# Removed materials import - now using user-defined crystal structures

app = Flask(__name__)
app.config['SECRET_KEY'] = 'displath-gui-secret'
socketio = SocketIO(app, cors_allowed_origins="*")

# Store running processes
running_processes = {}
output_queues = {}

def read_output(process, process_id, output_queue):
    """Read output from Julia process and put it in queue"""
    while True:
        line = process.stdout.readline()
        if not line and process.poll() is not None:
            break
        if line:
            output_queue.put(line.strip())
            socketio.emit('output', {'process_id': process_id, 'data': line.strip()}, room=process_id)
    
    # Process finished
    running_processes.pop(process_id, None)
    output_queues.pop(process_id, None)
    socketio.emit('process_finished', {'process_id': process_id}, room=process_id)

@app.route('/')
def index():
    return render_template('index.html')

def generate_julia_script(params):
    """Generate Julia script from parameters"""
    # Build type dictionary for crystal elements
    type_dict_lines = []
    next_type_id = 1
    
    # Parse crystal elements
    crystal_elements = params.get('crystal_elements', [])
    for elem in crystal_elements:
        if elem.get('symbol'):
            ed = float(elem.get("ed", 20.0))
            displacement_energy = float(elem.get("displacement_energy", 10.0))
            type_dict_lines.append(f'    {next_type_id} => Element("{elem["symbol"]}", {ed}, {displacement_energy})')
            next_type_id += 1
    
    # Add ion type
    ion_symbol = params.get('ion_symbol', 'Ne')
    ion_displacement = float(params.get('ion_displacement_energy', 0.1))
    ion_type_id = next_type_id
    type_dict_lines.append(f'    {ion_type_id} => Element("{ion_symbol}", 0.1, {ion_displacement})')
    
    # Get parameters
    ion_energy = params.get('ion_energy', 1000)
    repetitions = params.get('repetitions', 1000)
    irradiation_mode = params.get('irradiation_mode', 'independent')
    
    # Dynamic or static load
    is_dynamic = params.get('dynamic_load', False)
    is_amorphous = params.get('amorphous_mode', False)
    
    script = f'''home = ENV["ARCS_HOME"]
const IS_DYNAMIC_LOAD = {str(is_dynamic).lower()}
include(home * "/src/DISPLATH.jl")

# Materials & box
primaryVectors = [{params.get("primary_vectors", "1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0")}]
boxSizes = [{params.get("box_sizes", "10.0, 10.0, 10.0")}]
periodic = [{params.get("periodic", "true, true, true")}]
latticeRanges = [{params.get("lattice_ranges", "0 10; 0 10; 0 10")}]
basis = [{params.get("basis", "0.0 0.0 0.0")}]
basisTypes = [{params.get("basis_types", "1")}]
inputGridVectors = [{params.get("grid_vectors", "2.0 0.0 0.0; 0.0 2.0 0.0; 0.0 0.0 2.0")}]

typeDict = Dict(
{(',' + chr(10)).join(type_dict_lines)}
)

# Parameters
pMax = {float(params.get("pmax", 1.45))}
vacancyRecoverDistance = {float(params.get("vacancy_recover_distance", 1.0))}
seed = {int(params.get("seed", 42))}; const THREAD_RNG = [StableRNG(seed + t) for t in 1:Threads.nthreads()]
stopEnergy = {float(params.get("stop_energy", 0.1))}
'''

    # Always include temperature parameters
    script += f'''temperature = {float(params.get("temperature", 0.0))}
DebyeTemperature = {float(params.get("debye_temperature", 645.0))}

# Parameters constructor always includes temperature and DebyeTemperature
'''
    
    if is_amorphous:
        script += f'''parameters = Parameters(primaryVectors, latticeRanges, basisTypes, basis, pMax, vacancyRecoverDistance, typeDict;
                        stopEnergy=stopEnergy, temperature=temperature, DebyeTemperature=DebyeTemperature, amorphous=true)
'''
    else:
        script += f'''parameters = Parameters(primaryVectors, latticeRanges, basisTypes, basis, pMax, vacancyRecoverDistance, typeDict;
                        stopEnergy=stopEnergy, temperature=temperature, DebyeTemperature=DebyeTemperature)
'''

    script += f'''
# Process
simulator = Simulator(boxSizes, inputGridVectors, parameters)
@dump "init.dump" simulator.atoms

function Irradiation(simulator::Simulator, energy::Float64, iteration::Int)
    {f'Restore!(simulator)' if irradiation_mode == 'independent' else '# Continuous mode - no restore'}
    # Ion position: random x,y within box range, z at box_max_z - 0.1
    box_x_range = latticeRanges[1,2] - latticeRanges[1,1]
    box_y_range = latticeRanges[2,2] - latticeRanges[2,1]
    box_max_z = latticeRanges[3,2]
    
    ionPosition = [rand() * box_x_range + latticeRanges[1,1], 
                   rand() * box_y_range + latticeRanges[2,1], 
                   box_max_z + 0.1]
    ion = Atom({ion_type_id}, ionPosition, parameters)
    SetVelocityDirection!(ion, [{float(params.get("ion_direction_x", 0.0))}, {float(params.get("ion_direction_y", 0.0))}, {float(params.get("ion_direction_z", -1.0))}])
    SetEnergy!(ion, energy)
    push!(simulator, ion)
    
    # Store ion index for R_p calculation
    ion_index = length(simulator.atoms)
    
    Cascade!(ion, simulator)
'''
    
    # Add vacancy recording if enabled
    if params.get('record_vacancy', True):
        script += f'''    
    # Calculate vacancy count
    _, Vs = DefectStatics(simulator)
    nV = length(Vs)
    @record "nV.csv" "$(iteration),{ion_energy},$(nV)" "iteration,energy,nV"
'''
    
    # Add R_p recording if enabled
    if params.get('record_penetration', True):
        script += f'''    
    # Calculate penetration depth R_p
    ion_final_z = simulator.atoms[ion_index].coordinate[3]
    R_p = box_max_z - ion_final_z
    @record "R_p.csv" "$(iteration),{ion_energy},$(R_p)" "iteration,energy,R_p"
'''
    
    script += f'''    
    return nothing
end

{'Save!(simulator)' if irradiation_mode == 'independent' else '# Continuous mode - initial state saved once'}
@showprogress for i in 1:{repetitions}
    nV = Irradiation(simulator, {float(ion_energy)}, i)
end
@dump "final.dump" simulator.atoms
'''
    
    return script

@app.route('/run_simulation', methods=['POST'])
def run_simulation():
    """Run a Julia simulation"""
    data = request.json
    
    # Check if we need to generate script from parameters
    if 'parameters' in data:
        julia_script = generate_julia_script(data['parameters'])
    else:
        julia_script = data.get('script', '')
    
    num_threads = data.get('threads', 1)
    
    # Generate unique process ID
    process_id = f"process_{int(time.time() * 1000)}"
    
    # Get working directory from parameters, default to current directory
    working_dir = data.get('parameters', {}).get('working_directory', os.getcwd()) if 'parameters' in data else os.getcwd()
    
    # Ensure working directory exists
    try:
        os.makedirs(working_dir, exist_ok=True)
    except Exception as e:
        return jsonify({
            'success': False,
            'error': f'Cannot create working directory: {str(e)}'
        })
    
    # Create Julia file in working directory
    julia_file = os.path.join(working_dir, f"main_{process_id}.jl")
    with open(julia_file, 'w') as f:
        f.write(julia_script)
    
    try:
        # Start Julia process in working directory
        cmd = ['julia', f'-t{num_threads}', julia_file]
        process = subprocess.Popen(
            cmd,
            stdout=subprocess.PIPE,
            stderr=subprocess.STDOUT,
            universal_newlines=True,
            bufsize=1,
            cwd=working_dir  # Set working directory for Julia process
        )
        
        # Store process and create output queue
        running_processes[process_id] = process
        output_queues[process_id] = queue.Queue()
        
        # Start thread to read output
        thread = threading.Thread(
            target=read_output,
            args=(process, process_id, output_queues[process_id])
        )
        thread.daemon = True
        thread.start()
        
        return jsonify({
            'success': True,
            'process_id': process_id,
            'script': julia_script  # Return generated script for preview
        })
    
    except Exception as e:
        return jsonify({
            'success': False,
            'error': str(e)
        })

@app.route('/stop_simulation/<process_id>', methods=['POST'])
def stop_simulation(process_id):
    """Stop a running simulation"""
    if process_id in running_processes:
        process = running_processes[process_id]
        process.terminate()
        return jsonify({'success': True})
    return jsonify({'success': False, 'error': 'Process not found'})

@app.route('/get_element_presets', methods=['GET'])
def get_element_presets():
    """Get common element presets for quick input"""
    element_presets = {
        'C': {'symbol': 'C', 'ed': 22.0, 'displacement_energy': 11.0},
        'Si': {'symbol': 'Si', 'ed': 20.0, 'displacement_energy': 10.0},
        'B': {'symbol': 'B', 'ed': 18.3, 'displacement_energy': 9.15},
        'N': {'symbol': 'N', 'ed': 18.6, 'displacement_energy': 9.3},
        'S': {'symbol': 'S', 'ed': 7.8, 'displacement_energy': 3.9},
        'Mo': {'symbol': 'Mo', 'ed': 29.1, 'displacement_energy': 14.55},
        'He': {'symbol': 'He', 'ed': 0.1, 'displacement_energy': 0.1},
        'Ne': {'symbol': 'Ne', 'ed': 0.1, 'displacement_energy': 0.1},
        'Ar': {'symbol': 'Ar', 'ed': 0.1, 'displacement_energy': 0.1},
        'Kr': {'symbol': 'Kr', 'ed': 0.1, 'displacement_energy': 0.1},
        'Xe': {'symbol': 'Xe', 'ed': 0.1, 'displacement_energy': 0.1}
    }
    return jsonify(element_presets)

@app.route('/get_examples', methods=['GET'])
def get_examples():
    """Get list of example scripts"""
    examples_dir = os.path.join(os.path.dirname(__file__), '../../..', 'examples')
    examples = []
    
    for root, dirs, files in os.walk(examples_dir):
        for file in files:
            if file == 'main.jl':
                rel_path = os.path.relpath(os.path.join(root, file), examples_dir)
                examples.append(rel_path)
    
    return jsonify(examples)

@app.route('/load_example/<path:example_path>', methods=['GET'])
def load_example(example_path):
    """Load an example script"""
    examples_dir = os.path.join(os.path.dirname(__file__), '../../..', 'examples')
    file_path = os.path.join(examples_dir, example_path)
    
    try:
        with open(file_path, 'r') as f:
            content = f.read()
        return jsonify({
            'success': True,
            'content': content
        })
    except Exception as e:
        return jsonify({
            'success': False,
            'error': str(e)
        })

@socketio.on('connect')
def handle_connect():
    print('Client connected')

@socketio.on('disconnect')
def handle_disconnect():
    print('Client disconnected')

@app.route('/get_current_directory', methods=['GET'])
def get_current_directory():
    """Get current working directory"""
    try:
        current_dir = os.getcwd()
        return jsonify({
            'success': True,
            'directory': current_dir
        })
    except Exception as e:
        return jsonify({
            'success': False,
            'error': str(e)
        })

@app.route('/create_directory', methods=['POST'])
def create_directory():
    """Create and validate directory"""
    data = request.json
    directory = data.get('directory', '').strip()
    
    if not directory:
        return jsonify({
            'success': False,
            'error': 'Directory path cannot be empty'
        })
    
    try:
        # Expand user path (~ to home directory)
        directory = os.path.expanduser(directory)
        
        # Create directory if it doesn't exist
        os.makedirs(directory, exist_ok=True)
        
        # Verify it's accessible
        if not os.path.isdir(directory):
            return jsonify({
                'success': False,
                'error': 'Failed to create directory'
            })
        
        if not os.access(directory, os.W_OK):
            return jsonify({
                'success': False,
                'error': 'Directory is not writable'
            })
        
        return jsonify({
            'success': True,
            'directory': directory
        })
    
    except Exception as e:
        return jsonify({
            'success': False,
            'error': str(e)
        })

@app.route('/list_directory', methods=['POST'])
def list_directory():
    """List directory contents for autocomplete"""
    data = request.json
    path = data.get('path', '').strip()
    
    try:
        # Handle empty path - start from home or current directory
        if not path:
            path = os.path.expanduser('~')
        else:
            # Expand user path
            path = os.path.expanduser(path)
        
        # Split into directory and partial name for filtering
        if os.path.isdir(path) and not path.endswith(os.sep):
            # If it's a directory and doesn't end with separator, add it
            directory = path
            partial = ''
        else:
            # Split path into directory and partial filename
            directory = os.path.dirname(path)
            partial = os.path.basename(path)
        
        # If directory is empty (e.g., when typing just "test"), use current directory
        if not directory:
            directory = os.getcwd()
        
        # List directory contents
        if not os.path.exists(directory) or not os.path.isdir(directory):
            return jsonify({
                'success': True,
                'suggestions': []
            })
        
        suggestions = []
        try:
            items = os.listdir(directory)
            for item in items:
                # Filter by partial match (case-insensitive)
                if partial and not item.lower().startswith(partial.lower()):
                    continue
                
                full_path = os.path.join(directory, item)
                is_dir = os.path.isdir(full_path)
                
                # Only include directories and skip hidden files unless explicitly typed
                if not partial.startswith('.') and item.startswith('.'):
                    continue
                
                if is_dir:
                    suggestions.append({
                        'name': item,
                        'path': full_path + os.sep,
                        'type': 'directory',
                        'display': item + '/'
                    })
        except PermissionError:
            # Can't read directory, return empty suggestions
            pass
        
        # Sort suggestions: directories first, then alphabetically
        suggestions.sort(key=lambda x: (x['type'] != 'directory', x['name'].lower()))
        
        # Limit suggestions to prevent UI overload
        suggestions = suggestions[:20]
        
        return jsonify({
            'success': True,
            'suggestions': suggestions,
            'current_directory': directory
        })
    
    except Exception as e:
        return jsonify({
            'success': False,
            'error': str(e),
            'suggestions': []
        })

@socketio.on('connect')
def handle_connect():
    print('Client connected')

@socketio.on('disconnect')
def handle_disconnect():
    print('Client disconnected')

@socketio.on('join_process')
def handle_join_process(data):
    """Join a process room to receive its output"""
    process_id = data.get('process_id')
    if process_id:
        join_room(process_id)

if __name__ == '__main__':
    socketio.run(app, debug=True, host='0.0.0.0', port=5001, allow_unsafe_werkzeug=True)