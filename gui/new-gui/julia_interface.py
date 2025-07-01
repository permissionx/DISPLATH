"""
Safe Julia interface for DISPLATH GUI
Handles Julia integration with proper security and error handling
"""

import os
import sys
import json
import tempfile
import subprocess
from pathlib import Path
from typing import Dict, List, Any, Optional, Tuple
import time
from logging_config import get_logger, log_simulation_event, log_performance_metric

logger = get_logger(__name__)

class JuliaError(Exception):
    """Custom exception for Julia-related errors"""
    pass

class JuliaInterface:
    """Safe interface to Julia DISPLATH code"""
    
    def __init__(self):
        self.julia_path = None
        self.displath_path = None
        self.initialized = False
        self._find_julia_paths()
        
    def _find_julia_paths(self):
        """Find Julia executable and DISPLATH source paths"""
        try:
            # Find Julia executable
            julia_cmd = subprocess.run(['which', 'julia'], 
                                     capture_output=True, text=True)
            if julia_cmd.returncode == 0:
                self.julia_path = julia_cmd.stdout.strip()
                logger.info(f"Found Julia at: {self.julia_path}")
            else:
                logger.error("Julia executable not found in PATH")
                return
            
            # Find DISPLATH source directory
            current_dir = Path(__file__).parent
            displath_src = current_dir.parent.parent / "src" / "main.jl"
            
            if displath_src.exists():
                self.displath_path = str(displath_src)
                logger.info(f"Found DISPLATH at: {self.displath_path}")
                self.initialized = True
            else:
                logger.error(f"DISPLATH main.jl not found at: {displath_src}")
                
        except Exception as e:
            logger.error(f"Error finding Julia paths: {e}")
            
    def validate_parameters(self, params: Dict[str, Any]) -> Tuple[bool, List[str]]:
        """Validate simulation parameters before running"""
        errors = []
        
        try:
            # Validate cell vectors
            if 'cellVectors' not in params:
                errors.append("Cell vectors are required")
            else:
                cv = params['cellVectors']
                if not isinstance(cv, list) or len(cv) != 3:
                    errors.append("Cell vectors must be 3x3 matrix")
                else:
                    for i, row in enumerate(cv):
                        if not isinstance(row, list) or len(row) != 3:
                            errors.append(f"Cell vector row {i+1} must have 3 elements")
                        else:
                            for j, val in enumerate(row):
                                if not isinstance(val, (int, float)) or val < 0:
                                    errors.append(f"Cell vector [{i+1},{j+1}] must be positive number")
            
            # Validate box ranges
            if 'boxRanges' not in params:
                errors.append("Box ranges are required")
            else:
                br = params['boxRanges']
                if not isinstance(br, list) or len(br) != 3:
                    errors.append("Box ranges must have 3 elements")
                else:
                    for i, val in enumerate(br):
                        if not isinstance(val, int) or val <= 0:
                            errors.append(f"Box range {i+1} must be positive integer")
            
            # Validate basis atoms
            if 'basisAtoms' not in params or len(params['basisAtoms']) == 0:
                errors.append("At least one basis atom is required")
            else:
                for i, atom in enumerate(params['basisAtoms']):
                    required_fields = ['x', 'y', 'z', 'type']
                    for field in required_fields:
                        if field not in atom:
                            errors.append(f"Basis atom {i+1} missing {field}")
                        elif not isinstance(atom[field], (int, float)):
                            errors.append(f"Basis atom {i+1} {field} must be number")
            
            # Validate lattice ranges
            if 'latticeRanges' in params:
                lr = params['latticeRanges']
                if isinstance(lr, list) and len(lr) == 3:
                    for i, range_dict in enumerate(lr):
                        if not isinstance(range_dict, dict):
                            errors.append(f"Lattice range {i+1} must be object with min/max")
                        elif 'min' not in range_dict or 'max' not in range_dict:
                            errors.append(f"Lattice range {i+1} must have min and max")
                        elif range_dict['min'] >= range_dict['max']:
                            errors.append(f"Lattice range {i+1} min must be less than max")
            
            # Validate numerical parameters
            numerical_params = {
                'temperature': (0, 10000, "Temperature must be between 0-10000 K"),
                'debyeTemperature': (0, 10000, "Debye temperature must be between 0-10000 K"),
                'stopEnergy': (0.1, 1000000, "Stop energy must be between 0.1-1000000 eV"),
                'pMax': (0.1, 100, "pMax must be between 0.1-100 Å"),
                'vacancyRecoverDistance': (0, 100, "Vacancy recover distance must be between 0-100 Å"),
                'ionEnergy': (1, 10000000, "Ion energy must be between 1-10000000 eV"),
                'numRuns': (1, 100000, "Number of runs must be between 1-100000")
            }
            
            for param, (min_val, max_val, error_msg) in numerical_params.items():
                if param in params:
                    val = params[param]
                    if not isinstance(val, (int, float)) or val < min_val or val > max_val:
                        errors.append(error_msg)
            
        except Exception as e:
            logger.error(f"Error validating parameters: {e}")
            errors.append(f"Parameter validation error: {str(e)}")
        
        return len(errors) == 0, errors
    
    def create_julia_script(self, params: Dict[str, Any], script_path: str):
        """Create safe Julia script with validated parameters"""
        
        # Sanitize and validate all parameters
        is_valid, errors = self.validate_parameters(params)
        if not is_valid:
            raise JuliaError(f"Invalid parameters: {'; '.join(errors)}")
        
        # Extract and sanitize parameters
        cell_vectors = params['cellVectors']
        box_ranges = params['boxRanges']
        basis_atoms = params['basisAtoms']
        lattice_ranges = params.get('latticeRanges', [
            {'min': 0, 'max': box_ranges[0]},
            {'min': 0, 'max': box_ranges[1]},
            {'min': 0, 'max': box_ranges[2]}
        ])
        
        # Build Julia script safely
        script_content = f'''
# DISPLATH Simulation Script - Auto-generated
using Pkg
using LinearAlgebra

# Include DISPLATH
include("{self.displath_path}")

function run_simulation()
    try
        println("Starting DISPLATH simulation...")
        
        # Primary vectors (unit cell)
        primaryVectors = [
            {cell_vectors[0][0]} {cell_vectors[0][1]} {cell_vectors[0][2]};
            {cell_vectors[1][0]} {cell_vectors[1][1]} {cell_vectors[1][2]};
            {cell_vectors[2][0]} {cell_vectors[2][1]} {cell_vectors[2][2]}
        ]
        
        # Lattice ranges
        latticeRanges = [
            {lattice_ranges[0]['min']} {lattice_ranges[0]['max']};
            {lattice_ranges[1]['min']} {lattice_ranges[1]['max']};
            {lattice_ranges[2]['min']} {lattice_ranges[2]['max']}
        ]
        
        # Basis types and positions
        basisTypes = [{', '.join(str(atom['type']) for atom in basis_atoms)}]
        basis = [
{chr(10).join(f"            {atom['x']} {atom['y']} {atom['z']};" for atom in basis_atoms)}
        ]
        
        # Simulation parameters
        θτRepository = "{params.get('thetaTauRepository', '../../thetatau_repository')}"
        pMax = {params.get('pMax', 1.45)}
        vacancyRecoverDistance = {params.get('vacancyRecoverDistance', 0.0)}
        temperature = {params.get('temperature', 300.0)}
        debyeTemperature = {params.get('debyeTemperature', 645.0)}
        stopEnergy = {params.get('stopEnergy', 20.0)}
        
        # Create type dictionary
        typeDict = Dict{{Int64, Element}}()
        
        # Default elements (Si, C, N)
        typeDict[1] = Element("Si", 1.46, 28.0855, 14, 22.2, 11.1, 0.0, 0.0)
        typeDict[2] = Element("C", 1.46, 12.0107, 6, 22.2, 11.1, 0.0, 0.0)
        typeDict[3] = Element("N", 1.46, 14.0067, 7, 22.2, 11.1, 0.0, 0.0)
        
        # Create parameters
        parameters = Parameters(
            primaryVectors, latticeRanges, basisTypes, basis,
            θτRepository, pMax, vacancyRecoverDistance, typeDict;
            temperature = temperature,
            DebyeTemperature = debyeTemperature,
            stopEnergy = stopEnergy,
            isDumpInCascade = {str(params.get('isDumpInCascade', False)).lower()},
            isLog = true
        )
        
        println("Parameters created successfully")
        
        # Create simulator
        boxSizes = [{', '.join(map(str, box_ranges))}]
        gridVectors = [6.0 0.0 0.0; 0.0 6.0 0.0; 0.0 0.0 6.0]
        
        simulator = Simulator(boxSizes, gridVectors, parameters)
        println("Simulator created successfully")
        
        # Save initial state
        Save!(simulator)
        println("Initial state saved")
        
        # Run cascades
        numRuns = {params.get('numRuns', 100)}
        totalVacancies = 0
        ionType = {params.get('ionType', 1)}
        ionEnergy = {params.get('ionEnergy', 100000.0)}
        
        ionPosition = [{', '.join(map(str, params.get('ionPosition', [852.0, 852.0, 2047.0])))}]
        ionDirection = [{', '.join(map(str, params.get('ionDirection', [0.0, 0.0, -1.0])))}]
        
        println("Starting $numRuns cascades...")
        
        for i in 1:numRuns
            # Create ion
            ion = Atom(ionType, ionPosition, parameters)
            SetEnergy!(ion, ionEnergy)
            SetVelocityDirection!(ion, ionDirection)
            push!(simulator, ion)
            
            println("Run $i: Ion created with energy $ionEnergy eV")
            
            # Run cascade
            Cascade!(ion, simulator)
            
            # Count vacancies
            vacancies = CountVacancies(simulator)
            totalVacancies += vacancies
            
            println("Run $i completed: $vacancies vacancies")
            
            # Restore for next run
            if i < numRuns
                Restore!(simulator)
            end
            
            # Progress update
            progress = round(i / numRuns * 100, digits=1)
            println("Progress: $progress%")
        end
        
        println("Simulation completed!")
        println("Total vacancies: $totalVacancies")
        println("Average vacancies per cascade: $(totalVacancies / numRuns)")
        
        # Write results to file
        results = Dict(
            "totalVacancies" => totalVacancies,
            "averageVacancies" => totalVacancies / numRuns,
            "numRuns" => numRuns,
            "ionEnergy" => ionEnergy,
            "completed" => true
        )
        
        open("simulation_results.json", "w") do f
            write(f, JSON.json(results))
        end
        
        return results
        
    catch e
        println("ERROR: $e")
        
        error_results = Dict(
            "error" => string(e),
            "completed" => false
        )
        
        open("simulation_results.json", "w") do f
            write(f, JSON.json(error_results))
        end
        
        rethrow(e)
    end
end

# Run the simulation
results = run_simulation()
println("Script completed")
'''
        
        # Write script to file
        with open(script_path, 'w') as f:
            f.write(script_content)
        
        logger.info(f"Julia script created: {script_path}")
    
    def run_simulation(self, params: Dict[str, Any], sim_id: str) -> Dict[str, Any]:
        """Run DISPLATH simulation safely"""
        
        if not self.initialized:
            raise JuliaError("Julia interface not properly initialized")
        
        start_time = time.time()
        
        try:
            log_simulation_event(sim_id, "start", "Starting simulation", 
                                {"parameters": params})
            
            # Create temporary directory for this simulation
            with tempfile.TemporaryDirectory() as temp_dir:
                script_path = os.path.join(temp_dir, f"sim_{sim_id}.jl")
                results_path = os.path.join(temp_dir, "simulation_results.json")
                
                # Create Julia script
                self.create_julia_script(params, script_path)
                
                # Run Julia script
                cmd = [self.julia_path, script_path]
                
                log_simulation_event(sim_id, "julia_start", 
                                   f"Executing Julia command: {' '.join(cmd)}")
                
                # Execute with timeout
                timeout = params.get('timeout', 3600)  # 1 hour default
                process = subprocess.run(
                    cmd,
                    cwd=temp_dir,
                    capture_output=True,
                    text=True,
                    timeout=timeout
                )
                
                # Log Julia output
                if process.stdout:
                    log_simulation_event(sim_id, "julia_stdout", 
                                       "Julia output", {"output": process.stdout})
                
                if process.stderr:
                    log_simulation_event(sim_id, "julia_stderr", 
                                       "Julia errors", {"errors": process.stderr})
                
                # Check for errors
                if process.returncode != 0:
                    error_msg = f"Julia process failed with code {process.returncode}"
                    if process.stderr:
                        error_msg += f": {process.stderr}"
                    raise JuliaError(error_msg)
                
                # Read results
                if os.path.exists(results_path):
                    with open(results_path, 'r') as f:
                        results = json.load(f)
                else:
                    raise JuliaError("No results file generated")
                
                # Log completion
                duration = time.time() - start_time
                log_performance_metric("simulation_run", duration, 
                                     {"simulation_id": sim_id, "results": results})
                
                log_simulation_event(sim_id, "complete", "Simulation completed successfully",
                                   {"results": results, "duration": duration})
                
                return results
                
        except subprocess.TimeoutExpired:
            error_msg = f"Simulation timed out after {timeout} seconds"
            log_simulation_event(sim_id, "timeout", error_msg)
            raise JuliaError(error_msg)
            
        except Exception as e:
            duration = time.time() - start_time
            log_simulation_event(sim_id, "error", f"Simulation failed: {str(e)}",
                               {"duration": duration})
            raise JuliaError(f"Simulation failed: {str(e)}")

# Global Julia interface instance
julia_interface = JuliaInterface()