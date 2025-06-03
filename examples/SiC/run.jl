home = "/beegfs/home/xuke/Researches/Irradiation_Li-Tianzhao/4.DISPLATH/DISPLATH/"
include(home * "src/main.jl")


# Box and atoms 
a = 4.36
primaryVectors = [a 0.0 0.0; 0.0 a 0.0; 0.0 0.0 a]
boxSizes = [50, 50, 1400]
inputGridVectors = [a*2.1 0.0 0.0; 0.0 a*2.1 0.0; 0.0 0.0 a*2.1]
latticeRanges = [0 50; 0 50; 2 1200]
basis = [0.0 0.0 0.0; 0.25 0.25 0.25]
basisTypes = [1, 2]


# Parameters
θτRepository = home * "thetatau_repository/"
pMax = 4.0
vacancyRecoverDistance = 4.0
typeDict = Dict(
    1 => Element("Si", 43.0, 22.0),  
    2 => Element("C", 40.0, 20.0),
    3 => Element("N", 1.0, 1.0)    
)
temperature = 300.0
DebyeTemperature = 490.0
parameters = Parameters(θτRepository, pMax, vacancyRecoverDistance, typeDict; 
                        temperature=temperature, DebyeTemperature=DebyeTemperature)

# Run
simulator = Simulator(primaryVectors, boxSizes, inputGridVectors, latticeRanges, basis, basisTypes, parameters)

Save!(simulator)
energy = 10000.0
@showprogress for _ in 1:100
    ionPosition = RandomPointInCircle(20.0) + [110.0, 110.0, 5240.0] 
    ion = Atom(3, ionPosition, parameters)
    SetVelocityDirection!(ion, [0.,0.,-1.])
    SetEnergy!(ion,energy)
    push!(simulator, ion)
    Cascade!(ion, simulator)    
end


vacancies, interstitials = DefectStatics(simulator)
Natoms = simulator.atoms[simulator.atomNumberWhenStore:end]
Nzs = [atom.coordinate[3] for atom in Natoms]
open("Nzs.csv", "w") do file
    write(file, "z\n")
    for z in Nzs
        write(file, "$z\n")
    end
end





