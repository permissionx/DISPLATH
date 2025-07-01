from flask import Flask, jsonify, request, send_from_directory
from flask_cors import CORS
from flask_limiter import Limiter
from flask_limiter.util import get_remote_address
import threading
import time
from datetime import datetime, timedelta
import uuid
import os
from werkzeug.middleware.proxy_fix import ProxyFix
from typing import Dict, Any, Optional, List
import weakref
import gc

# Import our secure modules
from logging_config import get_logger, log_api_request, log_simulation_event, log_performance_metric
from julia_interface import julia_interface, JuliaError

app = Flask(__name__)
app.wsgi_app = ProxyFix(app.wsgi_app, x_for=1, x_proto=1, x_host=1, x_prefix=1)

# Configure CORS properly
CORS(app, 
     origins=["http://localhost:3000", "http://localhost:5000", "http://127.0.0.1:5000"],
     methods=["GET", "POST", "PUT", "DELETE"],
     allow_headers=["Content-Type", "Authorization"])

# Rate limiting
limiter = Limiter(
    key_func=get_remote_address,
    default_limits=["100 per hour", "10 per minute"]
)
limiter.init_app(app)

# Logger
logger = get_logger(__name__)

# Thread-safe simulation management
simulations = {}  # type: Dict[str, SimulationRunner]
simulation_lock = threading.RLock()
cleanup_lock = threading.Lock()

# Background cleanup thread
cleanup_thread = None
shutdown_event = threading.Event()

class SimulationRunner:
    """Thread-safe simulation runner with proper resource management"""
    
    def __init__(self, sim_id: str, parameters: Dict[str, Any]):
        self.sim_id = sim_id
        self.parameters = parameters.copy()  # Defensive copy
        self.status = "initializing"
        self.progress = 0
        self.results = {}
        self.logs = []
        self.start_time = datetime.now()
        self.thread: Optional[threading.Thread] = None
        self.cancelled = threading.Event()
        self._lock = threading.Lock()
        
        # Validate parameters on creation
        is_valid, errors = julia_interface.validate_parameters(parameters)
        if not is_valid:
            raise ValueError(f"Invalid parameters: {'; '.join(errors)}")
        
        logger.info(f"Created simulation runner {sim_id}")
        
    def log(self, message: str, level: str = "info"):
        """Thread-safe logging"""
        with self._lock:
            timestamp = datetime.now().strftime("%H:%M:%S")
            log_entry = {
                "time": timestamp,
                "message": message,
                "level": level
            }
            self.logs.append(log_entry)
            
            # Keep only last 1000 log entries to prevent memory issues
            if len(self.logs) > 1000:
                self.logs = self.logs[-1000:]
        
        # Also log to main logger
        log_simulation_event(self.sim_id, "log", message, {"level": level})
        
    def get_status(self) -> Dict[str, Any]:
        """Get current status safely"""
        with self._lock:
            return {
                "simulationId": self.sim_id,
                "status": self.status,
                "progress": self.progress,
                "logs": self.logs[-10:],  # Last 10 entries
                "results": self.results.copy() if self.status == "completed" else None,
                "startTime": self.start_time.isoformat(),
                "duration": (datetime.now() - self.start_time).total_seconds()
            }
    
    def cancel(self):
        """Cancel the simulation"""
        with self._lock:
            if self.status in ["running", "initializing"]:
                self.cancelled.set()
                self.status = "cancelled"
                self.log("Simulation cancelled by user", "warning")
                logger.info(f"Simulation {self.sim_id} cancelled")
                
    def run(self):
        """Run the simulation safely"""
        try:
            if self.cancelled.is_set():
                return
                
            with self._lock:
                self.status = "running"
            
            self.log("Starting DISPLATH simulation...")
            
            # Run Julia simulation
            start_time = time.time()
            results = julia_interface.run_simulation(self.parameters, self.sim_id)
            duration = time.time() - start_time
            
            if self.cancelled.is_set():
                with self._lock:
                    self.status = "cancelled"
                return
                
            # Update results
            with self._lock:
                self.results = results
                self.status = "completed"
                
            self.log(f"Simulation completed successfully in {duration:.1f}s")
            log_performance_metric("full_simulation", duration, 
                                 {"simulation_id": self.sim_id})
            
        except JuliaError as e:
            with self._lock:
                self.status = "error"
                self.results = {"error": str(e)}
            self.log(f"Julia error: {str(e)}", "error")
            logger.error(f"Simulation {self.sim_id} failed with Julia error: {e}")
            
        except Exception as e:
            with self._lock:
                self.status = "error" 
                self.results = {"error": f"Unexpected error: {str(e)}"}
            self.log(f"Unexpected error: {str(e)}", "error")
            logger.error(f"Simulation {self.sim_id} failed with unexpected error: {e}", 
                        exc_info=True)
    
    def start_async(self):
        """Start simulation in background thread"""
        if self.thread is not None:
            raise RuntimeError("Simulation already started")
            
        self.thread = threading.Thread(target=self.run, name=f"sim-{self.sim_id}")
        self.thread.daemon = True
        self.thread.start()
        logger.info(f"Started background thread for simulation {self.sim_id}")

def cleanup_old_simulations():
    """Background task to cleanup old simulations"""
    while not shutdown_event.is_set():
        try:
            with cleanup_lock:
                current_time = datetime.now()
                cutoff_time = current_time - timedelta(hours=24)  # Keep for 24 hours
                
                with simulation_lock:
                    to_remove = []
                    for sim_id, runner in simulations.items():
                        # Remove old completed/failed simulations
                        if (runner.start_time < cutoff_time and 
                            runner.status in ["completed", "error", "cancelled"]):
                            to_remove.append(sim_id)
                    
                    for sim_id in to_remove:
                        del simulations[sim_id]
                        logger.info(f"Cleaned up old simulation {sim_id}")
                
                # Force garbage collection
                gc.collect()
                
        except Exception as e:
            logger.error(f"Error in cleanup task: {e}")
        
        # Wait 1 hour between cleanups
        shutdown_event.wait(3600)

# Error handlers
@app.errorhandler(404)
def not_found(error):
    log_api_request(request.endpoint or "unknown", request.method, user_ip=request.remote_addr)
    return jsonify({"error": "Not found"}), 404

@app.errorhandler(400)
def bad_request(error):
    log_api_request(request.endpoint or "unknown", request.method, user_ip=request.remote_addr)
    return jsonify({"error": "Bad request"}), 400

@app.errorhandler(500)
def internal_error(error):
    logger.error(f"Internal server error: {error}", exc_info=True)
    return jsonify({"error": "Internal server error"}), 500

@app.before_request
def before_request():
    """Log all API requests"""
    if request.path.startswith('/api/'):
        log_api_request(request.endpoint or request.path, request.method, 
                       request.get_json(silent=True), request.remote_addr)

@app.route('/')
def index():
    """Serve main page"""
    return send_from_directory('.', 'index.html')

@app.route('/<path:path>')
def serve_static(path):
    """Serve static files"""
    try:
        return send_from_directory('.', path)
    except FileNotFoundError:
        return jsonify({"error": "File not found"}), 404

@app.route('/api/health', methods=['GET'])
def health_check():
    """Health check endpoint"""
    return jsonify({
        "status": "healthy",
        "timestamp": datetime.now().isoformat(),
        "julia_available": julia_interface.initialized,
        "active_simulations": len(simulations)
    })

@app.route('/api/validate', methods=['POST'])
@limiter.limit("20 per minute")
def validate_parameters():
    """Validate simulation parameters"""
    try:
        params = request.get_json()
        if not params:
            return jsonify({"error": "No parameters provided"}), 400
        
        is_valid, errors = julia_interface.validate_parameters(params)
        
        return jsonify({
            "valid": is_valid,
            "errors": errors
        })
        
    except Exception as e:
        logger.error(f"Error validating parameters: {e}", exc_info=True)
        return jsonify({"error": "Validation failed"}), 500

@app.route('/api/simulate', methods=['POST'])
@limiter.limit("5 per minute")
def start_simulation():
    """Start a new simulation"""
    try:
        params = request.get_json()
        if not params:
            return jsonify({"error": "No parameters provided"}), 400
        
        # Generate unique simulation ID
        sim_id = str(uuid.uuid4())
        
        with simulation_lock:
            # Check if we have too many active simulations
            active_count = sum(1 for runner in simulations.values() 
                             if runner.status in ["initializing", "running"])
            
            if active_count >= 10:  # Max 10 concurrent simulations
                return jsonify({"error": "Too many active simulations"}), 429
            
            # Create and start simulation runner
            try:
                runner = SimulationRunner(sim_id, params)
                simulations[sim_id] = runner
                runner.start_async()
                
            except ValueError as e:
                return jsonify({"error": str(e)}), 400
            except Exception as e:
                logger.error(f"Error creating simulation: {e}")
                return jsonify({"error": "Failed to create simulation"}), 500
        
        logger.info(f"Started simulation {sim_id}")
        return jsonify({
            "simulationId": sim_id,
            "status": "started"
        })
        
    except Exception as e:
        logger.error(f"Error starting simulation: {e}", exc_info=True)
        return jsonify({"error": "Failed to start simulation"}), 500

@app.route('/api/simulations', methods=['GET'])
def list_simulations():
    """List all simulations"""
    try:
        with simulation_lock:
            sims = []
            for sim_id, runner in simulations.items():
                sims.append({
                    "simulationId": sim_id,
                    "status": runner.status,
                    "startTime": runner.start_time.isoformat(),
                    "progress": runner.progress
                })
        
        return jsonify({"simulations": sims})
        
    except Exception as e:
        logger.error(f"Error listing simulations: {e}")
        return jsonify({"error": "Failed to list simulations"}), 500

@app.route('/api/simulation/<sim_id>/status', methods=['GET'])
def get_simulation_status(sim_id):
    """Get simulation status and progress"""
    try:
        with simulation_lock:
            if sim_id not in simulations:
                return jsonify({"error": "Simulation not found"}), 404
            
            runner = simulations[sim_id]
            status = runner.get_status()
        
        return jsonify(status)
        
    except Exception as e:
        logger.error(f"Error getting simulation status: {e}")
        return jsonify({"error": "Failed to get status"}), 500

@app.route('/api/simulation/<sim_id>/logs', methods=['GET'])
def get_simulation_logs(sim_id):
    """Get full simulation logs"""
    try:
        with simulation_lock:
            if sim_id not in simulations:
                return jsonify({"error": "Simulation not found"}), 404
            
            runner = simulations[sim_id]
            with runner._lock:
                logs = runner.logs.copy()
        
        # Pagination
        page = request.args.get('page', 1, type=int)
        per_page = request.args.get('per_page', 100, type=int)
        
        start_idx = (page - 1) * per_page
        end_idx = start_idx + per_page
        
        return jsonify({
            "logs": logs[start_idx:end_idx],
            "total": len(logs),
            "page": page,
            "per_page": per_page
        })
        
    except Exception as e:
        logger.error(f"Error getting simulation logs: {e}")
        return jsonify({"error": "Failed to get logs"}), 500

@app.route('/api/simulation/<sim_id>/cancel', methods=['POST'])
def cancel_simulation(sim_id):
    """Cancel a running simulation"""
    try:
        with simulation_lock:
            if sim_id not in simulations:
                return jsonify({"error": "Simulation not found"}), 404
            
            runner = simulations[sim_id]
            runner.cancel()
        
        logger.info(f"Cancelled simulation {sim_id}")
        return jsonify({
            "simulationId": sim_id,
            "status": "cancelled"
        })
        
    except Exception as e:
        logger.error(f"Error cancelling simulation: {e}")
        return jsonify({"error": "Failed to cancel simulation"}), 500

@app.route('/api/simulation/<sim_id>', methods=['DELETE'])
def delete_simulation(sim_id):
    """Delete a simulation and its data"""
    try:
        with simulation_lock:
            if sim_id not in simulations:
                return jsonify({"error": "Simulation not found"}), 404
            
            runner = simulations[sim_id]
            
            # Can only delete completed, cancelled, or failed simulations
            if runner.status in ["running", "initializing"]:
                return jsonify({"error": "Cannot delete running simulation"}), 400
            
            del simulations[sim_id]
        
        logger.info(f"Deleted simulation {sim_id}")
        return jsonify({"message": "Simulation deleted"})
        
    except Exception as e:
        logger.error(f"Error deleting simulation: {e}")
        return jsonify({"error": "Failed to delete simulation"}), 500

@app.route('/api/materials/presets', methods=['GET'])
def get_material_presets():
    """Get predefined material configurations"""
    try:
        presets = {
            "silicon": {
                "name": "Silicon (Diamond)",
                "cellVectors": [[5.431, 0.0, 0.0], [0.0, 5.431, 0.0], [0.0, 0.0, 5.431]],
                "basisAtoms": [
                    {"x": 0.0, "y": 0.0, "z": 0.0, "type": 1},
                    {"x": 0.25, "y": 0.25, "z": 0.25, "type": 1},
                    {"x": 0.5, "y": 0.5, "z": 0.0, "type": 1},
                    {"x": 0.75, "y": 0.75, "z": 0.25, "type": 1},
                    {"x": 0.5, "y": 0.0, "z": 0.5, "type": 1},
                    {"x": 0.75, "y": 0.25, "z": 0.75, "type": 1},
                    {"x": 0.0, "y": 0.5, "z": 0.5, "type": 1},
                    {"x": 0.25, "y": 0.75, "z": 0.75, "type": 1}
                ],
                "boxRanges": [100, 100, 100],
                "latticeRanges": [
                    {"min": 0, "max": 100},
                    {"min": 0, "max": 100},
                    {"min": 0, "max": 100}
                ]
            },
            "graphene": {
                "name": "Graphene",
                "cellVectors": [[4.26, 0.0, 0.0], [0.0, 4.263, 0.0], [0.0, 0.0, 20.0]],
                "basisAtoms": [
                    {"x": 0.0, "y": 0.0, "z": 0.0, "type": 2},
                    {"x": 0.333, "y": 0.667, "z": 0.0, "type": 2}
                ],
                "boxRanges": [400, 400, 20],
                "latticeRanges": [
                    {"min": 0, "max": 400},
                    {"min": 0, "max": 400},
                    {"min": 2, "max": 3}
                ]
            },
            "mos2": {
                "name": "MoS₂ (2H)",
                "cellVectors": [[3.16, 0.0, 0.0], [0.0, 3.16, 0.0], [0.0, 0.0, 12.3]],
                "basisAtoms": [
                    {"x": 0.0, "y": 0.0, "z": 0.0, "type": 1},
                    {"x": 0.333, "y": 0.667, "z": 0.25, "type": 2},
                    {"x": 0.333, "y": 0.667, "z": -0.25, "type": 2}
                ],
                "boxRanges": [200, 200, 50],
                "latticeRanges": [
                    {"min": 0, "max": 200},
                    {"min": 0, "max": 200},
                    {"min": 0, "max": 50}
                ]
            }
        }
        
        return jsonify(presets)
        
    except Exception as e:
        logger.error(f"Error getting material presets: {e}")
        return jsonify({"error": "Failed to get presets"}), 500

def start_cleanup_thread():
    """Start background cleanup thread"""
    global cleanup_thread
    
    if cleanup_thread is None or not cleanup_thread.is_alive():
        cleanup_thread = threading.Thread(target=cleanup_old_simulations, 
                                         name="cleanup", daemon=True)
        cleanup_thread.start()
        logger.info("Started cleanup thread")

def shutdown_handler():
    """Cleanup on shutdown"""
    logger.info("Shutting down server...")
    shutdown_event.set()
    
    # Cancel all running simulations
    with simulation_lock:
        for runner in simulations.values():
            if runner.status in ["running", "initializing"]:
                runner.cancel()
    
    logger.info("Shutdown complete")

import atexit
atexit.register(shutdown_handler)

if __name__ == '__main__':
    try:
        logger.info("Starting DISPLATH GUI server...")
        
        # Start cleanup thread
        start_cleanup_thread()
        
        # Check Julia availability
        if not julia_interface.initialized:
            logger.warning("Julia interface not initialized - simulations will fail")
        else:
            logger.info("Julia interface ready")
        
        logger.info("Server starting on http://localhost:5000")
        app.run(debug=False, port=5000, host='0.0.0.0', threaded=True)
        
    except KeyboardInterrupt:
        logger.info("Server stopped by user")
    except Exception as e:
        logger.error(f"Server startup failed: {e}", exc_info=True)
    finally:
        shutdown_handler()