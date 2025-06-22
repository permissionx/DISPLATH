home = "/beegfs/home/xuke/Researches/Irradiation_Li-Tianzhao/4.DISPLATH/DISPLATH/"
include(home * "/src/main.jl")


using Random

# Box
a = 1.42
b = 3.35
primaryVectors = [3.0*a 0.0 0.0; 0.0 3.0^0.5*a 0.0; 0.0 0.0 b]
boxSizes = [10, 20, 10]
inputGridVectors = [a*2.1 0.0 0.0; 0.0 a*2.1 0.0; 0.0 0.0 a*2.1]  # never be same as primaryVectors 
latticeRanges = [0 10; 0 20; 5 6]  

basis = [0.0 0.0 0.0; 1.0/3.0 0.0 0.0; 1.0/2.0 1.0/2.0 0.0; 5.0/6.0 1.0/2.0 0.0]
basisTypes = [1, 1, 1, 1]

# Parameters
θτRepository = home * "thetatau_repository/"
pMax = 1.45
vacancyRecoverDistance = 1.3
typeDict = Dict(
    1 => Element("C", 22.0, 22.0),  
    2 => Element("Ne", 22.0, 22.0)  
)
seed = 43
const THREAD_RNG = [StableRNG(seed + t) for t in 1:Threads.nthreads()]
parameters = Parameters(primaryVectors, latticeRanges, basisTypes, basis,
                        θτRepository, pMax, vacancyRecoverDistance, typeDict, seed)


# Initialize
simulator = Simulator_dynamicLoad(boxSizes, inputGridVectors, parameters) 
Save!(simulator)



function Irradiation(simulator::Simulator, energy::Float64)
    #if i % 100 == 0
        #println("Irradiation time: ", i)
    #end
    Restore!(simulator)
    ionPosition = RandomInAnUnitGrapheneCell(1.42) + [18.855045, 22.981482313368623, 20]
    ion = Atom(2, ionPosition, parameters)
    SetVelocityDirection!(ion, [0.0, 0.0, -1.0])
    SetEnergy!(ion,energy)
    push!(simulator, ion)
    Cascade!(ion, simulator)
end


Random.seed!(42)
computerNumberPerEnergy = 10000
vacancy_data = Dict{Int64, Vector{Int64}}()
energy = 500.0
n = 0
@showprogress  for i in 1:computerNumberPerEnergy
    Irradiation(simulator, energy)
end

