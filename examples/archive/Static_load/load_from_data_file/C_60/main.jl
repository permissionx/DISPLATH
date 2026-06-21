IS_DYNAMIC_LOAD = false   
home = ENV["ARCS_HOME"]
include(home * "/src/DISPLATH.jl")  

seed = 43   
const THREAD_RNG = [StableRNG(seed + t) for t in 1:Threads.nthreads()]

typeDict = Dict(               
    1 => Element("C", 22.0, 7.9),    #离位阈能， 结合能   单位：eV 
    2 => Element("Ne", 0.1, 0.1)     # 离子 0.1 
)

pMax = 3.5 
vacancyRecoverDistance = 0.0 

isDumpInCascade = false   # 输出辐照, default: false 
isNonQnl = true 
parameters = Parameters(pMax, vacancyRecoverDistance, typeDict;
                        isDumpInCascade = isDumpInCascade, isNonQnl = isNonQnl)

inputGridVectors = [4.0 0.0 0.0; 0.0 4.0 0.0; 0.0 0.0 4.0]   # cell 

fileName = "C60.data"
simulator = Simulator(fileName, inputGridVectors, parameters)
Save!(simulator)
@dump "C60.dump" simulator.atoms  #输出原子结构

for i in 1:100
    #Restore!(simulator) # 恢复到初始状态
    ionPosition = RandomPointInCircle(5.0)  + [10.2306,10.4537,15.8964]
    # 0,0,0 
    # ionPosition = [10.0, 10.0, 15.0]  # 离子位置
    ion = Atom(2, ionPosition, parameters)  
    SetVelocityDirection!(ion, [0.0, 0.0, -1.0])
    SetEnergy!(ion,1000.0) # 离子能量 eV 

    push!(simulator, ion)   
    Cascade!(ion, simulator)
    @dump "C60_after.dump" simulator.atoms
end







