home = "/beegfs/home/xuke/Researches/Irradiation_Li-Tianzhao/4.DISPLATH/DISPLATH/"
const BAlpha = 1.5
const IS_DYNAMIC_LOAD = true
include(home * "src/DISPLATH.jl")


# Box and atoms 
a = 5.431
primaryVectors = [a 0.0 0.0; 0.0 a 0.0; 0.0 0.0 a]
boxSizes = [400, 400, 1005]
inputGridVectors = [a*1.4 0.0 0.0; 0.0 a*1.4 0.0; 0.0 0.0 a*1.4]
latticeRanges = [0 400; 0 400; 2 1000]
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
pMax =  a/2/π^0.5
vacancyRecoverDistance = 0.0
typeDict = Dict(
    1 => Element("Si", 20.0, 10.0),  # Si element
    2 => Element("B", 1.0, 0.5)     # B for ion bombardment    
)
seed = 43
const THREAD_RNG = [StableRNG(seed + t) for t in 1:Threads.nthreads()]

temperature = 300.0
DebyeTemperature = 645.0
isDumpInCascade = false
stopEnergy = 20.0
nCascadeEveryLoad = 10
parameters = Parameters(primaryVectors, latticeRanges, basisTypes, basis,
                        θτRepository, pMax, vacancyRecoverDistance, typeDict; 
                        temperature=temperature, DebyeTemperature=DebyeTemperature, 
                        isDumpInCascade=isDumpInCascade, stopEnergy=stopEnergy, 
                        nCascadeEveryLoad=nCascadeEveryLoad)

# Run
simulator = Simulator(boxSizes, inputGridVectors, parameters)


for n in 1:10
    energy = 100_000.0 - (n-1) * 10_000.0
    for i in 1:1000
        Restore!(simulator)
        offset = [boxSizes[1]/2*a, boxSizes[2]/2*a, boxSizes[3]*a-2]
        rp = RandomPointInCircle(30.0)
        ionPosition = Vector{Float64}(rp) + offset
        ion = Atom(2, ionPosition, parameters)
        SetVelocityDirection!(ion, [0.12,0.,-1.])
        SetEnergy!(ion,energy)
        push!(simulator, ion)
        @time Cascade!(ion, simulator)
        z = simulator.atoms[1].coordinate[3]
        @record  "R_p.$(n).csv" "$(a*latticeRanges[3,2]-z)" isSmall=true
        @dump  "defects.$(n).dump" [simulator.atoms; simulator.vacancies] ["vx", "vy", "vz", "e"]
    end
end