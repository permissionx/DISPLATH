# Simple DISPLAΘ GUI launcher
using Pkg
Pkg.activate("..")
include("server.jl")
start_server(8080)