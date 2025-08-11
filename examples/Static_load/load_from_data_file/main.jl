IS_DYNAMIC_LOAD = false 
home = ENV["ARCS_HOME"]
include(home * "/src/DISPLATH.jl")

typeDict = Dict(
    1 => Element("C", 22.0, 11.0),  
    2 => Element("Ne", 0.1, 0.1)
)

pMax = 3.5 
vacancyRecoverDistance = 0.0 
parameters = Parameters(pMax, vacancyRecoverDistance, typeDict)

inputGridVectors = [1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0]
simulator = Simulator("C60.data", inputGridVectors, parameters)

@dump "C60.dump" simulator.atoms






