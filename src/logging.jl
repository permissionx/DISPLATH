# ===== Professional Logging Utilities =====

using Dates

# ANSI color codes
const GREEN = "\033[32m"
const YELLOW = "\033[33m"
const RED = "\033[31m"
const BLUE = "\033[34m"
const CYAN = "\033[36m"
const MAGENTA = "\033[35m"
const GRAY = "\033[90m"
const WHITE = "\033[37m"
const RESET = "\033[0m"
const BOLD = "\033[1m"
const DIM = "\033[2m"

"""
    log_info(args...)

Professional logging with timestamp and standard formatting.
"""
function log_info(args...)
    timestamp = Dates.format(now(), "HH:MM:SS.sss")
    print(DIM, GRAY, "[", timestamp, "] ", RESET)
    
    # Convert all arguments to strings and join them
    message = join(string.(args), "")
    print(message)
    println()
end

"""
    log_success(args...)

Success logging with timestamp and green checkmark.
"""
function log_success(args...)
    timestamp = Dates.format(now(), "HH:MM:SS.sss")
    print(DIM, GRAY, "[", timestamp, "] ", RESET)
    print(GREEN, "✓ ", RESET)
    
    message = join(string.(args), "")
    print(message)
    println()
end

"""
    log_warning(args...)

Warning logging with timestamp and yellow warning symbol.
"""
function log_warning(args...)
    timestamp = Dates.format(now(), "HH:MM:SS.sss")
    print(DIM, GRAY, "[", timestamp, "] ", RESET)
    print(YELLOW, "⚠ ", RESET)
    
    message = join(string.(args), "")
    print(YELLOW, message, RESET)
    println()
end

"""
    log_error(args...)

Error logging with timestamp and red X.
"""
function log_error(args...)
    timestamp = Dates.format(now(), "HH:MM:SS.sss")
    print(DIM, GRAY, "[", timestamp, "] ", RESET)
    print(RED, "✗ ", RESET)
    
    message = join(string.(args), "")
    print(RED, message, RESET)
    println()
end

"""
    log_debug(args...)

Debug logging with timestamp and dim gray color.
"""
function log_debug(args...)
    timestamp = Dates.format(now(), "HH:MM:SS.sss")
    print(DIM, GRAY, "[", timestamp, "] ", RESET)
    print(DIM, "• ", RESET)
    
    message = join(string.(args), "")
    print(DIM, message, RESET)
    println()
end

"""
    log_section(title)

Creates a section header with consistent formatting.
"""
function log_section(title::String)
    timestamp = Dates.format(now(), "HH:MM:SS.sss")
    print(DIM, GRAY, "[", timestamp, "] ", RESET)
    println()
    println(BOLD, "═══ ", title, " ═══", RESET)
end

"""
    log_separator()

Creates a simple separator line.
"""
function log_separator()
    println(DIM, "─" ^ 60, RESET)
end

# Export the logging functions
export log_info, log_success, log_warning, log_error, log_debug, log_section, log_separator 