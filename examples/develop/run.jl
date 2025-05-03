include("../../src/main.jl")
include("modules.jl")
using Random

a = 1.45
b = 3.0 
# Box and atoms
primaryVectors = [3.0*a 0.0 0.0; 0.0 3.0^0.5*a 0.0; 0.0 0.0 b]
boxSizes = [10,20,10]
inputGridVectors = [a*2.1 0.0 0.0; 0.0 a*2.1 0.0; 0.0 0.0 a*2.1]  # never be same as primaryVectors 
latticeRanges = [0 10; 0 20; 2 3]   
basis = [0.0 0.0 0.0; 1.0/3.0 0.0 0.0; 1.0/2.0 1.0/2.0 0.0; 5.0/6.0 1.0/2.0 0.0]
basisTypes = [1, 2, 1, 2]

# Parameters
θτRepository = "../../thetatau_repository/"
pMax = 1.45
vacancyRecoverDistance = 1.3
DTEMode = 2
DTEFile = "hBN.dte"
typeDict = Dict(
    1 => Element("B", 19.96, 19.96),  
    2 => Element("N", 22.77, 22.77),
    3 => Element("He", 1.0, 1.0)
)
parameters = Parameters(θτRepository, pMax, vacancyRecoverDistance, typeDict, DTEMode=DTEMode, DTEFile=DTEFile)


# Init simulator 
simulator = Simulator(primaryVectors, boxSizes, inputGridVectors, latticeRanges, basis, basisTypes, parameters)
Dump(simulator, "hBNH.dump", 0, false)
       
Random.seed!(40)
@showprogress for _ in 1:10000
    energy = 9000000.0
    simulator.nIrradiation += 1
    ionPosition = RandomInSquare(43.5, 50.2294) + [0.0, 0.0, 12.0]
    ion = Atom(3, ionPosition, parameters)
    SetVelocityDirection!(ion, [0.,0.,-1.])
    SetEnergy!(ion,energy) # Must be float number
    push!(simulator, ion)
    Cascade!(ion, simulator) 
    #if ion.isAlive
    #    delete!(simulator, ion)
    #end
    if simulator.nIrradiation % 100 == 0  
        Dump(simulator, "hBNH.dump", simulator.nIrradiation, true)  
    end
end

