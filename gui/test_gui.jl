# Test script for DISPLAΘ GUI
# This script will test the GUI step by step

println("=== DISPLAΘ GUI Test Script ===\n")

# Step 1: Check Julia environment
println("Step 1: Checking Julia environment...")
println("Julia version: ", VERSION)
println("Julia executable: ", Base.julia_cmd())
println("Working directory: ", pwd())
println()

# Step 2: Check required packages
println("Step 2: Checking required packages...")

required_packages = ["Gtk", "JSON3", "Dates", "LinearAlgebra", "QuadGK", "ProgressMeter", "StableRNGs", "Interpolations"]
missing_packages = String[]

for pkg in required_packages
    try
        eval(Meta.parse("using $pkg"))
        println("✓ $pkg - OK")
    catch e
        println("✗ $pkg - MISSING")
        push!(missing_packages, pkg)
    end
end

if !isempty(missing_packages)
    println("\n❌ Missing packages detected!")
    println("Run this command to install them:")
    println("julia -e 'using Pkg; Pkg.add([\"$(join(missing_packages, "\", \""))\"])'")
    println("\nOr use the install script:")
    println("julia gui/install_gui.jl")
    exit(1)
else
    println("\n✅ All required packages are available!")
end

# Step 3: Check main simulation code
println("\nStep 3: Testing main simulation code...")
try
    include("../src/main.jl")
    println("✓ Main simulation code loaded successfully")
    
    # Test basic types
    test_element = Element("Test", 1.0, 1.0)
    println("✓ Element creation works")
    
    # Test basic geometry
    test_vectors = [1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0]
    test_box = Box(test_vectors)
    println("✓ Box creation works")
    
catch e
    println("✗ Error loading main simulation code:")
    println(e)
    exit(1)
end

# Step 4: Test GUI components (without showing)
println("\nStep 4: Testing GUI components...")
try
    using Gtk
    
    # Test basic Gtk functionality
    window = GtkWindow("Test", 100, 100)
    println("✓ Basic window creation works")
    
    # Test widgets
    button = GtkButton("Test")
    entry = GtkEntry()
    label = GtkLabel("Test")
    println("✓ Basic widgets work")
    
    # Don't show the window, just test creation
    destroy(window)
    println("✓ GUI components test passed")
    
catch e
    println("✗ Error testing GUI components:")
    println(e)
    exit(1)
end

# Step 5: Test GUI file exists and syntax
println("\nStep 5: Testing GUI file...")
gui_file = "simple_gui.jl"
if isfile(gui_file)
    println("✓ GUI file exists: $gui_file")
    
    # Test syntax by parsing (not executing)
    try
        code = read(gui_file, String)
        parsed = Meta.parse("begin\n$code\nend")
        println("✓ GUI file syntax is valid")
    catch e
        println("✗ GUI file syntax error:")
        println(e)
        exit(1)
    end
else
    println("✗ GUI file not found: $gui_file")
    println("Make sure you're in the gui/ directory")
    exit(1)
end

# Step 6: Final recommendations
println("\n🎉 All tests passed!")
println("\n=== How to run the GUI ===")
println("1. Install missing packages (if any):")
println("   julia gui/install_gui.jl")
println()
println("2. Run the GUI:")
println("   julia gui/simple_gui.jl")
println()
println("3. Or from Julia REPL:")
println("   include(\"gui/simple_gui.jl\")")
println("   main()")
println()

# Step 7: Test with minimal parameters
println("=== Quick Simulation Test ===")
println("Testing with minimal parameters...")

try
    # Create minimal test parameters
    a = 4.36
    primaryVectors = [a 0.0 0.0; 0.0 a 0.0; 0.0 0.0 a]
    boxSizes = [5, 5, 10]  # Very small for testing
    inputGridVectors = [a*2.1 0.0 0.0; 0.0 a*2.1 0.0; 0.0 0.0 a*2.1]
    latticeRanges = [0 5; 0 5; 2 8]
    basis = [0.0 0.0 0.0; 0.25 0.25 0.25]
    basisTypes = [1, 2]
    
    θτRepository = "../thetatau_repository/"
    if !isdir(θτRepository)
        println("⚠️  Warning: θτRepository not found at $θτRepository")
        println("   GUI will work but simulations may fail")
    else
        println("✓ θτRepository found")
    end
    
    typeDict = Dict(
        1 => Element("Si", 43.0, 22.0),
        2 => Element("C", 40.0, 20.0),
        3 => Element("N", 1.0, 1.0)
    )
    
    parameters = Parameters(
        primaryVectors, latticeRanges, basisTypes, basis,
        θτRepository, 4.0, 4.0, typeDict, 12345;
        temperature=300.0, DebyeTemperature=490.0
    )
    
    println("✓ Parameter creation successful")
    
    # Test simulator creation (this might take a moment)
    println("Testing simulator creation (this may take a few seconds)...")
    simulator = Simulator(primaryVectors, boxSizes, inputGridVectors, latticeRanges, basis, basisTypes, parameters)
    println("✓ Simulator creation successful")
    
    println("✅ Core simulation components are working!")
    
catch e
    println("⚠️  Warning: Simulation test failed:")
    println(e)
    println("GUI interface will still work, but simulations may have issues")
end

println("\n🚀 Ready to test the GUI!")
println("Run: julia gui/simple_gui.jl")