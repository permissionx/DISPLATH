home = "/beegfs/home/xuke/Researches/Irradiation_Li-Tianzhao/4.DISPLATH/DISPLATH/"
include(home * "/src/main.jl")
include("modules.jl")
using Random

ncompute = 1
a = 1.42
b = 6.70
nlayers = [1,2,3,4,5]
nlayer = nlayers[ncompute]
# Box and atoms
primaryVectors = [3.0*a 0.0 0.0; 0.0 3.0^0.5*a 0.0; 0.0 0.0 b]
boxSizes = [10,20,10]
inputGridVectors = [a*2.2 0.0 0.0; 0.0 a*2.2 0.0; 0.0 0.0 a*2.2]  # never be same as primaryVectors 
latticeRanges = [0 10; 0 20; 2 2+nlayer]   # modify: layer number 
basis = [0.0 0.0 0.0; 1.0/3.0 0.0 0.0; 1.0/2.0 1.0/2.0 0.0; 5.0/6.0 1.0/2.0 0.0; 
         1.0/6.0 1.0/2.0 0.5; 1.0/3.0+1.0/6.0 0.0+1.0/2.0 0.5; 1.0/2.0+1.0/6.0 1.0/2.0+1.0/2.0 0.5; 5.0/6.0+1.0/6.0 1.0/2.0+1.0/2.0 0.5;]
basisTypes = [1, 1, 2, 1, 1, 2, 1, 1]

# Parameters
θτRepository = home * "/thetatau_repository/"
pMax = 1.42
vacancyRecoverDistance = 4.0
typeDict = Dict(
    1 => Element("C", 19.96, 19.96),
    2 => Element("C", 19.96, 19.96),
    3 => Element("Ne", 1.0, 1.0)      # modify: ion type  dte: 1.0 eV
)
seed = 42
const THREAD_RNG = [StableRNG(seed + t) for t in 1:Threads.nthreads()]
parameters = Parameters(primaryVectors, latticeRanges, basisTypes, basis,
                        θτRepository, pMax, vacancyRecoverDistance, typeDict, seed)
# Init simulator 
simulator = Simulator(boxSizes, inputGridVectors, parameters)
Save!(simulator)
@dump "graphite.dump" simulator.atoms
# exit()  


x1 = [10*1.2^x for x in 0:35]
x2 = [x1[end]*1.7^x for x in 1:15]
energys = [x1..., x2...]

vacancyDate = Dict{Int64, Tuple{Float64, Float64}}()


for (na, energy) in enumerate(energys)
    nVs = Vector{Int64}()
    nCount = 10000 # modify: optional  
    @showprogress desc="$(na)" for irun in 1:nCount
        nV = Irradiation(simulator, energy)
        push!(nVs, nV)
    end
    meanNVs = sum(nVs) / nCount
    vacancyDate[na] = (energy, meanNVs)
    simulator.nCascade = 0
end

open("vacancy.csv", "w") do f
    write(f, "n,energy,nVacancies\n")
    for n in sort(collect(keys(vacancyDate)))
        write(f, "$(n),$(vacancyDate[n][1]),$(vacancyDate[n][2])\n")
    end
end


