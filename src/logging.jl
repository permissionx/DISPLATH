# ===== Professional Logging Utilities =====

using Dates

# ANSI color codes
const GREEN = "\033[32m"
const YELLOW = "\033[33m"
const RED = "\033[31m"
const BLUE = "\033[34m"
const GRAY = "\033[90m"
const RESET = "\033[0m"
const BOLD = "\033[1m"
const DIM = "\033[2m"

"""
    log_info(args...)

Professional logging with timestamp and gray color formatting.
Replaces standard println with timestamped, colored output.
"""
function log_info(args...)
    timestamp = Dates.format(now(), "HH:MM:SS.sss")
    print(GRAY, "[", timestamp, "] ", RESET)
    
    # Convert all arguments to strings and join them
    message = join(string.(args), "")
    
    # Apply gray color to the message content
    print(GRAY, message, RESET)
    println()
end

"""
    log_success(args...)

Success logging with timestamp and green color.
"""
function log_success(args...)
    timestamp = Dates.format(now(), "HH:MM:SS.sss")
    print("\033[32m", "[", timestamp, "] ✅ ", RESET)
    
    message = join(string.(args), "")
    print("\033[32m", message, RESET)
    println()
end

"""
    log_warning(args...)

Warning logging with timestamp and yellow color.
"""
function log_warning(args...)
    timestamp = Dates.format(now(), "HH:MM:SS.sss")
    print("\033[33m", "[", timestamp, "] ⚠️ ", RESET)
    
    message = join(string.(args), "")
    print("\033[33m", message, RESET)
    println()
end

"""
    log_error(args...)

Error logging with timestamp and red color.
"""
function log_error(args...)
    timestamp = Dates.format(now(), "HH:MM:SS.sss")
    print("\033[31m", "[", timestamp, "] ❌ ", RESET)
    
    message = join(string.(args), "")
    print("\033[31m", message, RESET)
    println()
end

"""
    log_debug(args...)

Debug logging with timestamp and dim gray color.
"""
function log_debug(args...)
    timestamp = Dates.format(now(), "HH:MM:SS.sss")
    print(DIM, GRAY, "[", timestamp, "] 🐛 ", RESET)
    
    message = join(string.(args), "")
    print(DIM, GRAY, message, RESET)
    println()
end

# Export the logging functions
export log_info, log_success, log_warning, log_error, log_debug 