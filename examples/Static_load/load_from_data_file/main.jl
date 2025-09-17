IS_DYNAMIC_LOAD = false 
home = ENV["ARCS_HOME"]
include(home * "/src/DISPLATH.jl")

typeDict = Dict(
    1 => Element("C", 22.0, 11.0),  
    2 => Element("Ne", 0.1, 0.1)
)

pMax = 3.5 
vacancyRecoverDistance = 0.0 

isDumpInCascade = true
parameters = Parameters(pMax, vacancyRecoverDistance, typeDict;
                        isDumpInCascade = isDumpInCascade)

inputGridVectors = [4.0 0.0 0.0; 0.0 4.0 0.0; 0.0 0.0 4.0]

fileName = "C60.data"
simulator = Simulator(fileName, inputGridVectors, parameters)

@dump "C60.dump" simulator.atoms


ionPosition = [10.0, 10.0, 15.0]
ion = Atom(2, ionPosition, parameters)
SetVelocityDirection!(ion, [0.0, 0.0, -1.0])
SetEnergy!(ion,1000.0)
push!(simulator, ion)

Cascade!(ion, simulator)

@dump "C60_after.dump" simulator.atoms




