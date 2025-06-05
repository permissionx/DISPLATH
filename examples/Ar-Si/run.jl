home = "/beegfs/home/xuke/Researches/Irradiation_Li-Tianzhao/4.DISPLATH/DISPLATH/"
include(home * "src/main.jl")
using Random
using Profile

# Box and atoms 
a = 5.431
primaryVectors = [a 0.0 0.0; 0.0 a 0.0; 0.0 0.0 a]
boxSizes = [400, 400, 405]
inputGridVectors = [a*2.1 0.0 0.0; 0.0 a*2.1 0.0; 0.0 0.0 a*2.1]
latticeRanges = [0 400; 0 400; 0 400]
# Si diamond structure - conventional cell with 8 atoms
basis = [0.0 0.0 0.0;           # atom 1
         0.5 0.5 0.0;           # atom 2  
         0.5 0.0 0.5;           # atom 3
         0.0 0.5 0.5;           # atom 4
         0.25 0.25 0.25;        # atom 5
         0.75 0.75 0.25;        # atom 6
         0.75 0.25 0.75;        # atom 7
         0.25 0.75 0.75]        # atom 8
basisTypes = [1, 1, 1, 1, 1, 1, 1, 1]  # All 8 positions are Si atoms


# Parameters
θτRepository = home * "thetatau_repository/"
pMax = 4.0
vacancyRecoverDistance = 4.0
typeDict = Dict(
    1 => Element("Si", 22.5, 11.0),  # Si element
    2 => Element("Ar", 1.0, 1.0)     # B for ion bombardment    
)

temperature = 300.0
DebyeTemperature = 645.0
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


Random.seed!(43)
energy = 40000.0

println("Running 1 ion...")
if simulator.parameters.isDumpInCascade
    Dump(simulator, "Cascade.dump", 0, false)
end

ionPosition = RandomPointInCircle(30.0) + [1000.0, 1000.0, 2198.0] 
ion = Atom(2, ionPosition, parameters)
SetVelocityDirection!(ion, [0.,0.,-1.])
SetEnergy!(ion,energy)
push!(simulator, ion)
Log("E_atom,E_electron\n",type="w")
@time Cascade_dynamicLoad!(ion, simulator)
println("Initersitial number: $(length(simulator.atoms)); Vacancy number: $(length(simulator.vacancies))")








