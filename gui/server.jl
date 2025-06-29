using HTTP
using JSON3
using Dates
using Random

# Include the DISPLAΘ simulation code
home = dirname(dirname(@__FILE__))
const IS_DYNAMIC_LOAD = false # Will be set by GUI
include(joinpath(home, "src", "DISPLATH.jl"))

mutable struct SimulationConfig
    # Material parameters
    primaryVectors::Matrix{Float64}
    latticeRanges::Matrix{Int64}
    basisTypes::Vector{Int64}
    basis::Matrix{Float64}
    
    # Simulation parameters
    boxSizes::Vector{Int64}
    inputGridVectors::Matrix{Float64}
    pMax::Float64
    vacancyRecoverDistance::Float64
    temperature::Float64
    DebyeTemperature::Float64
    stopEnergy::Float64
    
    # Ion parameters
    ionType::Int64
    ionEnergy::Float64
    ionDirection::Vector{Float64}
    ionPosition::Vector{Float64}
    
    # Simulation settings
    nRuns::Int64
    isDynamicLoad::Bool
    nCascadeEveryLoad::Int64
    
    # Output settings
    outputName::String
    isDumpInCascade::Bool
    
    # Type dictionary
    typeDict::Dict{Int64, Element}
end

function create_default_config()
    # Default graphene configuration
    a = 1.42
    b = 6.70
    primaryVectors = [3.0*a 0.0 0.0; 0.0 3.0^0.5*a 0.0; 0.0 0.0 b]
    latticeRanges = [0 10; 0 20; 2 3]
    basis = [0.0 0.0 0.0; 1.0/3.0 0.0 0.0; 1.0/2.0 1.0/2.0 0.0; 5.0/6.0 1.0/2.0 0.0]
    basisTypes = [1, 1, 1, 1]
    
    boxSizes = [10, 20, 10]
    inputGridVectors = [a*2.2 0.0 0.0; 0.0 a*2.2 0.0; 0.0 0.0 a*2.2]
    
    typeDict = Dict(
        1 => Element("C", 19.96, 19.96),
        2 => Element("Ne", 1.0, 1.0)
    )
    
    return SimulationConfig(
        primaryVectors, latticeRanges, basisTypes, basis,
        boxSizes, inputGridVectors, 1.42, 4.0, 300.0, 1000.0, 10.0,
        2, 1000.0, [0.0, 0.0, -1.0], [0.0, 0.0, 0.0],
        100, false, 1,
        "simulation_output", false,
        typeDict
    )
end

function run_simulation(config::SimulationConfig)
    try
        # Set global dynamic load flag
        global IS_DYNAMIC_LOAD = config.isDynamicLoad
        
        # Create parameters
        θτRepository = joinpath(home, "thetatau_repository")
        seed = rand(1:10000)
        
        parameters = Parameters(
            config.primaryVectors, config.latticeRanges, config.basisTypes, config.basis,
            θτRepository, config.pMax, config.vacancyRecoverDistance, config.typeDict;
            temperature=config.temperature, DebyeTemperature=config.DebyeTemperature,
            stopEnergy=config.stopEnergy, isDumpInCascade=config.isDumpInCascade,
            nCascadeEveryLoad=config.nCascadeEveryLoad
        )
        
        # Initialize simulator
        simulator = Simulator(Vector{Int64}(config.boxSizes), config.inputGridVectors, parameters)
        Save!(simulator)
        
        # Run simulations
        vacancies = Int64[]
        
        for i in 1:config.nRuns
            Restore!(simulator)
            
            # Create ion
            ion = Atom(config.ionType, config.ionPosition, parameters)
            SetVelocityDirection!(ion, config.ionDirection)
            SetEnergy!(ion, config.ionEnergy)
            push!(simulator, ion)
            
            # Run cascade
            Cascade!(ion, simulator)
            
            # Count vacancies
            nV = CountVacancies(simulator)
            push!(vacancies, nV)
        end
        
        # Calculate statistics
        meanVacancies = sum(vacancies) / length(vacancies)
        stdVacancies = sqrt(sum((v - meanVacancies)^2 for v in vacancies) / length(vacancies))
        
        return Dict(
            "success" => true,
            "results" => Dict(
                "totalRuns" => config.nRuns,
                "meanVacancies" => meanVacancies,
                "stdVacancies" => stdVacancies,
                "allVacancies" => vacancies,
                "simulationTime" => now()
            )
        )
        
    catch e
        return Dict(
            "success" => false,
            "error" => string(e),
            "stacktrace" => stacktrace()
        )
    end
end

function handle_request(req::HTTP.Request)
    # CORS headers
    headers = [
        "Access-Control-Allow-Origin" => "*",
        "Access-Control-Allow-Methods" => "GET, POST, OPTIONS",
        "Access-Control-Allow-Headers" => "Content-Type",
        "Content-Type" => "application/json"
    ]
    
    if req.method == "OPTIONS"
        return HTTP.Response(200, headers, "")
    end
    
    try
        if req.target == "/" || req.target == "/index.html"
            # Serve the main HTML page
            html_content = read(joinpath(@__DIR__, "index.html"), String)
            return HTTP.Response(200, ["Content-Type" => "text/html"], html_content)
        elseif req.target == "/style.css"
            # Serve CSS
            css_content = read(joinpath(@__DIR__, "style.css"), String)
            return HTTP.Response(200, ["Content-Type" => "text/css"], css_content)
        elseif req.target == "/script.js"
            # Serve JavaScript
            js_content = read(joinpath(@__DIR__, "script.js"), String)
            return HTTP.Response(200, ["Content-Type" => "application/javascript"], js_content)
        elseif req.target == "/api/simulate" && req.method == "POST"
            # Handle simulation request
            data = JSON3.read(String(req.body))
            
            # Convert JSON data to SimulationConfig
            config = SimulationConfig(
                reshape(collect(data.primaryVectors), 3, 3)',
                reshape(collect(data.latticeRanges), 3, 2)',
                collect(data.basisTypes),
                reshape(collect(data.basis), length(data.basisTypes), 3)',
                collect(data.boxSizes),
                reshape(collect(data.inputGridVectors), 3, 3)',
                data.pMax, data.vacancyRecoverDistance,
                data.temperature, data.DebyeTemperature, data.stopEnergy,
                data.ionType, data.ionEnergy,
                collect(data.ionDirection), collect(data.ionPosition),
                data.nRuns, data.isDynamicLoad, data.nCascadeEveryLoad,
                data.outputName, data.isDumpInCascade,
                Dict(Int64(k) => Element(v.name, v.dte, v.bde) for (k,v) in data.typeDict)
            )
            
            result = run_simulation(config)
            return HTTP.Response(200, headers, JSON3.write(result))
            
        elseif req.target == "/api/presets"
            # Return preset configurations
            presets = Dict(
                "graphene" => Dict(
                    "name" => "Graphene Monolayer",
                    "primaryVectors" => [4.26, 0.0, 0.0, 0.0, 4.263, 0.0, 0.0, 0.0, 6.70],
                    "latticeRanges" => [0, 10, 0, 20, 2, 3],
                    "basis" => [0.0, 0.0, 0.0, 1.0/3.0, 0.0, 0.0, 1.0/2.0, 1.0/2.0, 0.0, 5.0/6.0, 1.0/2.0, 0.0],
                    "basisTypes" => [1, 1, 1, 1],
                    "typeDict" => Dict("1" => Dict("name" => "C", "dte" => 19.96, "bde" => 19.96), "2" => Dict("name" => "Ne", "dte" => 1.0, "bde" => 1.0))
                ),
                "silicon" => Dict(
                    "name" => "Silicon Crystal",
                    "primaryVectors" => [5.431, 0.0, 0.0, 0.0, 5.431, 0.0, 0.0, 0.0, 5.431],
                    "latticeRanges" => [0, 50, 0, 50, 2, 50],
                    "basis" => [0.0, 0.0, 0.0, 0.5, 0.5, 0.0, 0.5, 0.0, 0.5, 0.0, 0.5, 0.5, 0.25, 0.25, 0.25, 0.75, 0.75, 0.25, 0.75, 0.25, 0.75, 0.25, 0.75, 0.75],
                    "basisTypes" => [1, 1, 1, 1, 1, 1, 1, 1],
                    "typeDict" => Dict("1" => Dict("name" => "Si", "dte" => 20.0, "bde" => 10.0), "2" => Dict("name" => "B", "dte" => 1.0, "bde" => 0.5))
                )
            )
            return HTTP.Response(200, headers, JSON3.write(presets))
        else
            return HTTP.Response(404, headers, JSON3.write(Dict("error" => "Not found")))
        end
    catch e
        error_response = Dict("error" => string(e), "stacktrace" => string(stacktrace()))
        return HTTP.Response(500, headers, JSON3.write(error_response))
    end
end

function start_server(port::Int = 8080)
    println("Starting DISPLAΘ GUI server on port $port...")
    println("Open http://localhost:$port in your browser")
    
    server = HTTP.serve(handle_request, "0.0.0.0", port)
    return server
end

# Start the server if this script is run directly
if abspath(PROGRAM_FILE) == @__FILE__
    start_server()
end