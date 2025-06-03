home = "/beegfs/home/xuke/Researches/Irradiation_Li-Tianzhao/4.DISPLATH/DISPLATH/"
include(home * "src/main.jl")
using Random
using Profile

# Box and atoms 
a = 4.36
primaryVectors = [a 0.0 0.0; 0.0 a 0.0; 0.0 0.0 a]
boxSizes = [800, 800, 2002]
inputGridVectors = [a*4.1 0.0 0.0; 0.0 a*4.1 0.0; 0.0 0.0 a*4.1]
latticeRanges = [0 800; 0 800; 2 2000]
basis = [0.0 0.0 0.0; 0.25 0.25 0.25]
basisTypes = [1, 2]


# Parameters
θτRepository = home * "thetatau_repository/"
pMax = 4.0
vacancyRecoverDistance = 10.0
typeDict = Dict(
    1 => Element("Si", 20.0, 10.0),  
    2 => Element("C", 35.0, 17.5),
    3 => Element("N", 1.0, 1.0)    
)
temperature = 300.0
DebyeTemperature = 300.0
isDumpInCascade = true
isDynamicLoad = true
stopEnergy = 20.0
parameters = Parameters(primaryVectors, latticeRanges, basisTypes, basis,
                        θτRepository, pMax, vacancyRecoverDistance, typeDict; 
                        temperature=temperature, DebyeTemperature=DebyeTemperature, 
                        isDumpInCascade=isDumpInCascade, isDynamicLoad=isDynamicLoad,
                        stopEnergy=stopEnergy)

# Run
simulator = Simulator_dynamicLoad(boxSizes, inputGridVectors, parameters)
#cell = simulator.cellGrid.cells[10, 10, 10]
#@time LoadCell(cell, simulator)

Random.seed!(43)
energy = 28000.0
#ionPosition = RandomPointInCircle(20.0) + [110.0, 110.0, 5240.0] 
#ion = Atom(3, ionPosition, parameters)
#SetVelocityDirection!(ion, [0.,0.,-1.])
#SetEnergy!(ion,energy)
#push!(simulator, ion)
#Cascade_dynamicLoad!(ion, simulator)

#Profile.clear()
println("Running 10 ions...")
if simulator.parameters.isDumpInCascade
    Dump_dynamicLoad(simulator, "Cascade.dump", 0, false)
end
@showprogress for _ in 1:1000
    ionPosition = RandomPointInCircle(400.0) + [1740.0, 1740.0, 8726.0] 
    ion = Atom(3, ionPosition, parameters)
    SetVelocityDirection!(ion, [0.0,0.,-1.])
    SetEnergy!(ion,energy)
    push!(simulator, ion)
    Cascade_dynamicLoad!(ion, simulator)
end
#Profile.print(maxdepth=12, mincount=100)

#println("\n保存详细Profile结果到文件...")
#open("cascade_profile.txt", "w") do io
#    println(io, "=== Cascade_dynamicLoad Profile结果 ===")
#    println(io, "")
#    Profile.print(io, maxdepth=20, mincount=50)
#end







