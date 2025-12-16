IS_DYNAMIC_LOAD = false   
home = ENV["ARCS_HOME"]
include(home * "/src/DISPLATH.jl")  
seed = 43;const THREAD_RNG = [StableRNG(seed + t) for t in 1:Threads.nthreads()]

# Parameters
pMax = 3.5 
vacancyRecoverDistance = 0.0 
parameters = Parameters(pMax, vacancyRecoverDistance;
                        isDumpInCascade = true)

# Material
fileName = "C60.data"
typeDict = Dict(               
    1 => Element("C", 22.0, 7.9),    #离位阈能， 结合能   单位：eV 
    2 => Element("Ne", 0.1, 0.1)     # 离子 0.1 
)
inputGridVectors = [4.0 0.0 0.0; 0.0 4.0 0.0; 0.0 0.0 4.0]   # cell 
material = Material(fileName, typeDict, inputGridVectors, parameters;
                    replicate=[1,1,1])


simulator = Simulator(material, parameters)
Save!(simulator)
@dump "C60.dump" simulator.atoms  

ionPosition = RandomPointInCircle(5.0)  + [10.2306,10.4537,15.8964]
energy = 100.0
Irradiation!(simulator, energy, ionPosition, [0.0,0.0,-1.0], 2, parameters)
@dump "final.dump" simulator.atoms 

#=
function Irradiation(simulator::Simulator, energy::Float64)
    Restore!(simulator)
    ionPosition =  [12.0, 5.0, 17.0]
    ion = Atom(4, ionPosition, parameters)
    SetVelocityDirection!(ion, [0.0, 0.0, -1.0])
    SetEnergy!(ion,energy)
    push!(simulator, ion)
    Cascade!(ion, simulator)
end

Irradiation(simulator, 10.0)
=#