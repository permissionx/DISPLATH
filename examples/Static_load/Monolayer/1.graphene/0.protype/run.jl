# using BCA.jl
home = "/beegfs/home/xuke/Researches/Irradiation_Li-Tianzhao/4.DISPLATH/DISPLATH"
include(home * "/src/main.jl")
##include("modules.jl")


# Box
a = 1.42
b = 3.35
primaryVectors = [3.0*a 0.0 0.0; 0.0 3.0^0.5*a 0.0; 0.0 0.0 b]
boxSizes = [10, 20, 10]
inputGridVectors = [a*2.1 0.0 0.0; 0.0 a*2.1 0.0; 0.0 0.0 a*2.1]  # never be same as primaryVectors 
periodic = [true, true, false]

latticeRanges = [0 10; 0 20; 5 6]   
basis = [0.0 0.0 0.0; 1.0/3.0 0.0 0.0; 1.0/2.0 1.0/2.0 0.0; 5.0/6.0 1.0/2.0 0.0]
basisTypes = [1, 1, 1, 1]

# Parameters
θτRepository = "."
pMax = 1.45
vacancyRecoverDistance = 1.3
# for Ne ion 
typeDict = Dict(
    1 => Element("C", 22.0, 11.0),  
    2 => Element("Ne", 22.0, 11.0)  
)
seed = 42
EPowerRange = -1.0:0.045:8.0
pRange = 0.0:0.01:2.0
stopEnergy = 0.1
const THREAD_RNG = [StableRNG(seed + t) for t in 1:Threads.nthreads()]
parameters = Parameters(primaryVectors, latticeRanges, basisTypes, basis,
                        θτRepository, pMax, vacancyRecoverDistance, typeDict, seed;
                        EPowerRange=EPowerRange, pRange=pRange, stopEnergy=stopEnergy)


# Initialize
simulator = Simulator(boxSizes, inputGridVectors, parameters)  
Save!(simulator)

#Dump(simulator, "run.dump",0,false)


function Irradiation(simulator::Simulator, energy::Float64)
    #if i % 100 == 0
        #println("Irradiation time: ", i)
    #end
    Restore!(simulator)
    simulator.nCascade += 1
    ionPosition = RandomInAnUnitGrapheneCell(1.39667) + [18.855045, 22.981482313368623, 20]
    ion = Atom(2, ionPosition, parameters)
    SetVelocityDirection!(ion, [0.0, 0.0, -1.0])
    SetEnergy!(ion,energy)
    push!(simulator, ion)
    Cascade!(ion, simulator)
    Vs, _ = DefectStatics(simulator)
    return length(Vs)
end


Random.seed!(42)
computerNumberPerEnergy = 10000
vacancy_data = Dict{Int64, Vector{Int64}}()
x1 = [10*1.2^x for x in 0:35]
x2 = [x1[end]*1.7^x for x in 1:15]
energys = [x1..., x2...]
n = 0
for energy in energys
    global n += 1
    vacancyNumbers = Vector{Int64}(undef, computerNumberPerEnergy)
    @showprogress desc="In irradiation of energy order: $(n) ($(round(energy, digits=2)) eV)" for i in 1:computerNumberPerEnergy
        vacancyNumber = Irradiation(simulator, energy)
        vacancyNumbers[i] = vacancyNumber
    end
    vacancy_data[n] = vacancyNumbers
end


log_info("Outputting data...")
open("vacancy.csv", "w") do f
    write(f, "n,energy,nVacancies\n")
    for n in sort(collect(keys(vacancy_data)))
        for nVacancies in vacancy_data[n]
            write(f, "$(n),$(energys[n]),$(nVacancies)\n")
        end
    end
end

