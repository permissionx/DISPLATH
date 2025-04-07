include("../../src/main.jl")
a = 1.45
b = 3.0 

# Box and atoms
primaryVectors = [3.0*a 0.0 0.0; 0.0 3.0^0.5*a 0.0; 0.0 0.0 b]
boxSizes = [10,20,10]
inputGridVectors = [a*2.1 0.0 0.0; 0.0 a*2.1 0.0; 0.0 0.0 a*2.1]  # never be same as primaryVectors 
periodic = [true, true, false]
latticeRanges = [0 10; 0 20; 5 6]   
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
    3 => Element("Xe", 1.0, 1.0)
)
parameters = Parameters(θτRepository, pMax, vacancyRecoverDistance, typeDict, DTEMode=DTEMode, DTEFile=DTEFile)


# Init simulator 
simulator = Simulator(primaryVectors, boxSizes, inputGridVectors, periodic, latticeRanges, basis, basisTypes, parameters)


