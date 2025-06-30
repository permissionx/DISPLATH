#!/usr/bin/env julia

# Simple startup script for DISPLAΘ GUI
# Usage: julia start_gui.jl [port]

using Pkg
using Dates

# Activate the parent project environment
Pkg.activate("..")

# Include the server
include("server.jl")

# Include logging utilities
include("../src/logging.jl")

# Get port from command line argument or use default
port = length(ARGS) > 0 ? parse(Int, ARGS[1]) : 8080

log_info("=" ^ 60)
log_success("🚀 DISPLAΘ BCA Simulator Web GUI")
log_info("=" ^ 60)
log_info("")
log_info("📂 Project: $(basename(dirname(pwd())))")
log_info("🌐 Server:  http://localhost:$port")
log_info("⏰ Started: $(now())")
log_info("")
log_info("💡 Tips:")
log_info("   • Open the URL above in your web browser")
log_info("   • Use Ctrl+C to stop the server")
log_info("   • Check the console for error messages")
log_info("")
log_info("=" ^ 60)

try
    # Start the server
    start_server(port)
catch InterruptException
    log_info("\n👋 Server stopped by user")
catch e
    log_error("\n❌ Server error: $e")
    rethrow(e)
end