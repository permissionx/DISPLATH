# Installation script for DISPLAΘ GUI dependencies

println("Installing DISPLAΘ GUI dependencies...")

using Pkg

# Required packages for the GUI
required_packages = [
    "Gtk",           # GUI framework
    "JSON3",         # Parameter saving/loading
    "Dates"          # Timestamp functionality
]

println("Installing packages: $(join(required_packages, ", "))")

try
    # Install packages
    Pkg.add(required_packages)
    
    println("✓ All packages installed successfully!")
    println("\nTo run the GUI:")
    println("  julia gui/simple_gui.jl")
    println("\nOr from Julia REPL:")
    println("  include(\"gui/simple_gui.jl\")")
    println("  main()")
    
catch e
    println("✗ Error installing packages: $e")
    println("\nTry installing manually:")
    for pkg in required_packages
        println("  Pkg.add(\"$pkg\")")
    end
end