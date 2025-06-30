home = "/beegfs/home/xuke/Researches/Irradiation_Li-Tianzhao/4.DISPLATH/DISPLATH/"
const BAlpha = 1.5
const NAlpha = 1.5
const IS_DYNAMIC_LOAD = false   
include(home * "/src/DISPLATH.jl")

# Box
a = 1.45
b = 3.35
primaryVectors = [3.0*a 0.0 0.0; 0.0 3.0^0.5*a 0.0; 0.0 0.0 b]
basis = [0.0 0.0 0.0; 
         1.0/3.0 0.0 0.0; 
         1.0/2.0 1.0/2.0 0.0; 
         5.0/6.0 1.0/2.0 0.0]
basisTypes = [1, 2, 1, 2]
typeDict = Dict(
    1 => Element("B", 19.96, 19.96),  
    2 => Element("N", 22.77, 22.77),
    3 => Element("He", 0.1, 0.1)
)
lx = 50
ly = 30 
boxSizes = [lx, ly, 5]
inputGridVectors = [3.0 0.0 0.0; 0.0 3.0 0.0; 0.0 0.0 3.0]  # never be same as primaryVectors 
periodic = [true, true, false]
latticeRanges = [0 lx; 0 ly; 1 2]   
# Parameters
θτRepository = home * "thetatau_repository/"
pMax = 1.45
vacancyRecoverDistance = 1.3
seed = 42
const THREAD_RNG = [StableRNG(seed + t) for t in 1:Threads.nthreads()]
stopEnergy = 0.1
DTEMode = 2
DTEFile = "../../../dte_repository/hBN.dte"
parameters = Parameters(primaryVectors, latticeRanges, basisTypes, basis,
                        θτRepository, pMax, vacancyRecoverDistance, typeDict;
                        stopEnergy=stopEnergy, DTEMode=DTEMode, DTEFile=DTEFile, isDumpInCascade=false)
# Initialize
simulator = Simulator(boxSizes, inputGridVectors, parameters)  
@dump "hBN.dump" simulator.atoms  


function Irradiation(simulator::Simulator, energy::Float64)
    ionPosition = RandomPointInCircle(20.0) + [100.5037, 35.963, 10.0]
    ion = Atom(3, ionPosition, parameters)
    SetVelocityDirection!(ion, [0.0, 0.0, -1.0])
    SetEnergy!(ion,energy)
    push!(simulator, ion)
    Cascade!(ion, simulator)
end


computerNumberPerEnergy = 1  # custom
ens = [2.0]  #cutom 
Save!(simulator)
for en in ens
    energy = 10.0^en
    nVmean = 0 
    @showprogress desc="In irradiation of energy order: $(en) ($(round(energy, digits=2)) eV)" for i in 1:computerNumberPerEnergy
        Restore!(simulator)
        for j in 1:1000
            Irradiation(simulator, energy)
            if j % 10 == 0
                @dump "out/hBN.dump" simulator.atoms 
            end
        end
        Is, Vs = DefectStatics(simulator)
        nVmean += length(Vs) 
    end
    nVmean /= computerNumberPerEnergy
    @record "out/nV.csv" "$(energy),$(nVmean)"
end


