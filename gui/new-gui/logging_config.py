"""
Logging configuration for DISPLATH GUI
"""

import logging
import logging.handlers
import os
from datetime import datetime
import json

# Create logs directory if it doesn't exist
LOGS_DIR = os.path.join(os.path.dirname(__file__), 'logs')
os.makedirs(LOGS_DIR, exist_ok=True)

class JSONFormatter(logging.Formatter):
    """Custom JSON formatter for structured logging"""
    
    def format(self, record):
        log_data = {
            'timestamp': datetime.fromtimestamp(record.created).isoformat(),
            'level': record.levelname,
            'logger': record.name,
            'message': record.getMessage(),
            'module': record.module,
            'function': record.funcName,
            'line': record.lineno
        }
        
        # Add exception info if present
        if record.exc_info:
            log_data['exception'] = self.formatException(record.exc_info)
        
        # Add extra fields if present
        if hasattr(record, 'extra_fields'):
            log_data.update(record.extra_fields)
            
        return json.dumps(log_data, ensure_ascii=False)

def setup_logging():
    """Setup comprehensive logging configuration"""
    
    # Root logger
    root_logger = logging.getLogger()
    root_logger.setLevel(logging.DEBUG)
    
    # Clear any existing handlers
    root_logger.handlers.clear()
    
    # Console handler with color coding
    console_handler = logging.StreamHandler()
    console_handler.setLevel(logging.INFO)
    console_formatter = logging.Formatter(
        '%(asctime)s - %(name)s - %(levelname)s - %(message)s',
        datefmt='%Y-%m-%d %H:%M:%S'
    )
    console_handler.setFormatter(console_formatter)
    root_logger.addHandler(console_handler)
    
    # File handler for all logs (JSON format)
    all_logs_file = os.path.join(LOGS_DIR, 'displath_gui.log')
    file_handler = logging.handlers.RotatingFileHandler(
        all_logs_file,
        maxBytes=10*1024*1024,  # 10MB
        backupCount=5
    )
    file_handler.setLevel(logging.DEBUG)
    file_handler.setFormatter(JSONFormatter())
    root_logger.addHandler(file_handler)
    
    # Error file handler (JSON format)
    error_logs_file = os.path.join(LOGS_DIR, 'displath_errors.log')
    error_handler = logging.handlers.RotatingFileHandler(
        error_logs_file,
        maxBytes=5*1024*1024,  # 5MB
        backupCount=3
    )
    error_handler.setLevel(logging.ERROR)
    error_handler.setFormatter(JSONFormatter())
    root_logger.addHandler(error_handler)
    
    # Simulation logs file
    sim_logs_file = os.path.join(LOGS_DIR, 'simulations.log')
    sim_handler = logging.handlers.RotatingFileHandler(
        sim_logs_file,
        maxBytes=20*1024*1024,  # 20MB
        backupCount=10
    )
    sim_handler.setLevel(logging.INFO)
    sim_handler.setFormatter(JSONFormatter())
    
    # Create simulation logger
    sim_logger = logging.getLogger('simulation')
    sim_logger.addHandler(sim_handler)
    sim_logger.setLevel(logging.INFO)
    
    return root_logger

def get_logger(name):
    """Get a logger with the specified name"""
    return logging.getLogger(name)

def log_simulation_event(sim_id, event_type, message, extra_data=None):
    """Log simulation-specific events with structured data"""
    logger = logging.getLogger('simulation')
    
    extra_fields = {
        'simulation_id': sim_id,
        'event_type': event_type
    }
    
    if extra_data:
        extra_fields.update(extra_data)
    
    # Create a custom LogRecord with extra fields
    record = logger.makeRecord(
        logger.name, logging.INFO, __file__, 0, message, (), None
    )
    record.extra_fields = extra_fields
    
    logger.handle(record)

def log_api_request(endpoint, method, params=None, user_ip=None):
    """Log API requests for monitoring and debugging"""
    logger = logging.getLogger('api')
    
    extra_fields = {
        'endpoint': endpoint,
        'method': method,
        'user_ip': user_ip or 'unknown'
    }
    
    if params:
        # Don't log sensitive data
        safe_params = {k: v for k, v in params.items() 
                      if k not in ['password', 'token', 'secret']}
        extra_fields['parameters'] = safe_params
    
    record = logger.makeRecord(
        logger.name, logging.INFO, __file__, 0, 
        f"{method} {endpoint}", (), None
    )
    record.extra_fields = extra_fields
    
    logger.handle(record)

def log_performance_metric(operation, duration, extra_data=None):
    """Log performance metrics"""
    logger = logging.getLogger('performance')
    
    extra_fields = {
        'operation': operation,
        'duration_seconds': duration,
        'timestamp': datetime.now().isoformat()
    }
    
    if extra_data:
        extra_fields.update(extra_data)
    
    record = logger.makeRecord(
        logger.name, logging.INFO, __file__, 0,
        f"Performance: {operation} took {duration:.3f}s", (), None
    )
    record.extra_fields = extra_fields
    
    logger.handle(record)

# Initialize logging when module is imported
setup_logging()