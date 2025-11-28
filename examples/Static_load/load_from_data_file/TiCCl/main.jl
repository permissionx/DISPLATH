IS_DYNAMIC_LOAD = false   
home = ENV["ARCS_HOME"]
include(home * "/src/DISPLATH.jl")  

seed = 43   
const THREAD_RNG = [StableRNG(seed + t) for t in 1:Threads.nthreads()]

typeDict = Dict(               
    1 => Element("Ti", 22.0, 7.9),    #离位阈能， 结合能   单位：eV 
    2 => Element("C", 0.1, 0.1),     # 离子 0.1 
    3 => Element("Cl", 0.1, 0.1),     # 离子 0.1 
    4 => Element("He", 0.1, 0.1),     # 离子 0.1 
)

pMax = 3.5 
vacancyRecoverDistance = 0.0 
parameters = Parameters(pMax, vacancyRecoverDistance, typeDict;
                        isDumpInCascade = false, isNonQnl = true)

inputGridVectors = [4.0 0.0 0.0; 0.0 4.0 0.0; 0.0 0.0 4.0]   # cell 

fileName = "Ti2CCl.data"
simulator = Simulator(fileName, inputGridVectors, parameters; replicate=[5,5,1])
Save!(simulator)
@dump "TiCCl.dump" simulator.atoms  #输出原子结构








