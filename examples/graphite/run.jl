include("../../src/main.jl")
include("modules.jl")
using Random

a = 1.42
b = 6.70
# Box and atoms
primaryVectors = [3.0*a 0.0 0.0; 0.0 3.0^0.5*a 0.0; 0.0 0.0 b]
boxSizes = [10,20,10]
inputGridVectors = [a*2.1 0.0 0.0; 0.0 a*2.1 0.0; 0.0 0.0 a*2.1]  # never be same as primaryVectors 
latticeRanges = [0 10; 0 20; 2 5]   # modify: layer number 
basis = [0.0 0.0 0.0; 1.0/3.0 0.0 0.0; 1.0/2.0 1.0/2.0 0.0; 5.0/6.0 1.0/2.0 0.0; 
         1.0/6.0 1.0/2.0 0.5; 1.0/3.0+1.0/6.0 0.0+1.0/2.0 0.5; 1.0/2.0+1.0/6.0 1.0/2.0+1.0/2.0 0.5; 5.0/6.0+1.0/6.0 1.0/2.0+1.0/2.0 0.5;]
basisTypes = [1, 1, 2, 1, 1, 2, 1, 1]

# Parameters
θτRepository = "../../thetatau_repository/"
pMax = 1.42
vacancyRecoverDistance = 4.0
typeDict = Dict(
    1 => Element("C", 19.96, 19.96),
    2 => Element("C", 19.96, 19.96),
    3 => Element("Ne", 1.0, 1.0)      # modify: ion type  dte: 1.0 eV
)

parameters = Parameters(θτRepository, pMax, vacancyRecoverDistance, typeDict, DTEMode=4)
# Init simulator 
simulator = Simulator(primaryVectors, boxSizes, inputGridVectors, latticeRanges, basis, basisTypes, parameters)
Save!(simulator)
Dump(simulator, "graphite.dump", 0, false)
# exit()  

Random.seed!(42)  # modify: optional  
nVs = Vector{Int64}()


energy = 10000.0  # eV. Modify: ion energy 
nCount = 10000 # modify: optional  
@showprogress for irun in 1:nCount
    nV = Irradiation(simulator, energy)
    push!(nVs, nV)
    #if irun == 20
    #    Dump(simulator, "graphite_2.dump", 0, false)
    #end
end
meanNVs = sum(nVs) / nCount

open("Ne.$(energy).txt", "w") do file
    for nV in nVs
        write(file, "$(nV)\n")
    end
end

open("Ne.$(energy).mean.txt","w") do file
    write(file, "$(meanNVs)\n")
end

