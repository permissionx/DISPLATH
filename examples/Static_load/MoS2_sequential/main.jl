# using BCA.jl
home = "/beegfs/home/xuke/Researches/Irradiation_Li-Tianzhao/4.DISPLATH/DISPLATH/"
const IS_DYNAMIC_LOAD = false
include(home * "/src/DISPLATH.jl")
##include("modules.jl")


# Box
primaryVectors = [3.1962223053 0.0 0.0; 0.0 5.5360207558 0.0; 0.0 0.0 23.1298294067]
basis = [0.500000000 0.166667000 0.500000000;
         0.000000000 0.666667000 0.500000000;
         0.500000000 0.833333000 0.432137000;
         0.000000000 0.333333000 0.432137000;
         0.500000000 0.833333000 0.567863000;
         0.000000000 0.333333000 0.567863000;]
basisTypes = [1, 1, 2, 2, 2, 2]
typeDict = Dict(
    1 => Element("S", 7.8, 7.8/2),  
    2 => Element("Mo", 29.1, 29.1/2),
    3 => Element("Ar", 0.1, 0.1)
)
lx = 50
ly = 30 
boxSizes = [lx, ly, 3]
inputGridVectors = [6.0 0.0 0.0; 0.0 6.0 0.0; 0.0 0.0 6.0]  # never be same as primaryVectors 
periodic = [true, true, false]
latticeRanges = [0 lx; 0 ly; 1 2]   
# Parameters
θτRepository = home * "thetatau_repository/"

pMax = 1.45
vacancyRecoverDistance = 3.0
seed = 42
const THREAD_RNG = [StableRNG(seed + t) for t in 1:Threads.nthreads()]

EPowerRange = -1.0:0.045:8.0
pRange = 0.0:0.01:2.0
stopEnergy = 0.1
parameters = Parameters(primaryVectors, latticeRanges, basisTypes, basis,
                        θτRepository, pMax, vacancyRecoverDistance, typeDict;
                        EPowerRange=EPowerRange, pRange=pRange, stopEnergy=stopEnergy,
                        isDumpInCascade=false)


# Initialize
simulator = Simulator(boxSizes, inputGridVectors, parameters)  
@dump "MoS2.dump" simulator.atoms 


#
r = 67.0E-8 
S = pi * r^2
phi = 2.0E14     #  custom 
nIon = phi * S 
#

function Irradiation(simulator::Simulator, energy::Float64)
    ionPosition = RandomPointInCircle(67.0) + [81.5037, 83.963, 68.0]
    ion = Atom(3, ionPosition, parameters)
    SetVelocityDirection!(ion, [0.0, 0.0, -1.0])
    SetEnergy!(ion,energy)
    push!(simulator, ion)
    Cascade!(ion, simulator)
end


computerNumberPerEnergy = 1000  # custom

ens = 1:0.1:4  #cutom 

Save!(simulator)
for en in ens
    energy = 10^en
    nVmean = 0
    @showprogress desc="In irradiation of energy order: $(en) ($(round(energy, digits=2)) eV)" for i in 1:computerNumberPerEnergy
        Restore!(simulator)
        for j in 1:nIon
            Irradiation(simulator, energy)
        end
        Is, Vs = DefectStatics(simulator)
        nVmean += length(Vs)
    end
    nVmean /= computerNumberPerEnergy
    @record "nV.csv" "$(energy),$(nVmean)"
end



