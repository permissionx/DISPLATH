#!/usr/bin/env julia

# Simple startup script for DISPLAΘ GUI
# Usage: julia start_gui.jl [port]

using Pkg
using Dates

# Activate the parent project environment
Pkg.activate("..")

# Include the server
include("server.jl")

# Get port from command line argument or use default
port = length(ARGS) > 0 ? parse(Int, ARGS[1]) : 8080

println("=" ^ 60)
println("🚀 DISPLAΘ BCA Simulator Web GUI")
println("=" ^ 60)
println()
println("📂 Project: $(basename(dirname(pwd())))")
println("🌐 Server:  http://localhost:$port")
println("⏰ Started: $(now())")
println()
println("💡 Tips:")
println("   • Open the URL above in your web browser")
println("   • Use Ctrl+C to stop the server")
println("   • Check the console for error messages")
println()
println("=" ^ 60)

try
    # Start the server
    start_server(port)
catch InterruptException
    println("\n👋 Server stopped by user")
catch e
    println("\n❌ Server error: $e")
    rethrow(e)
end