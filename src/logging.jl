# ===== Professional Logging Utilities =====

using Dates

# ANSI color codes
const GRAY = "\033[90m"      # 灰色
const RESET = "\033[0m"      # 重置颜色
const BOLD = "\033[1m"       # 粗体
const DIM = "\033[2m"        # 暗淡

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
    print("\033[32m", "[", timestamp, "] ✅ ", RESET)  # 绿色
    
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
    print("\033[33m", "[", timestamp, "] ⚠️ ", RESET)  # 黄色
    
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
    print("\033[31m", "[", timestamp, "] ❌ ", RESET)  # 红色
    
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