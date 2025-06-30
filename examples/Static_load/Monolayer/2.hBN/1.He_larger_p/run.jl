# using BCA.jl
include("../../../../../src/main.jl")
include("modules.jl")


using Random


# Box and atoms 
a = 1.45
b = 3.35
primaryVectors = [3.0*a 0.0 0.0; 0.0 3.0^0.5*a 0.0; 0.0 0.0 b]
boxSizes = [20, 40, 10]
inputGridVectors = [a*3.1 0.0 0.0; 0.0 a*3.1 0.0; 0.0 0.0 a*3.1]  # never be same as primaryVectors 
latticeRanges = [0 20; 0 40; 2 3]   
basis = [0.0 0.0 0.0; 1.0/3.0 0.0 0.0; 1.0/2.0 1.0/2.0 0.0; 5.0/6.0 1.0/2.0 0.0]
basisTypes = [1, 2, 1, 2]

# Parameters
θτRepository = "../../../../../thetatau_repository/"
pMax = 3.1
vacancyRecoverDistance = 1.3
# for Xe ion 
typeDict = Dict(
    1 => Element("B", 19.96, 19.96),  
    2 => Element("N", 22.77, 22.77),
    3 => Element("He", 1.0, 1.0)
)


parameters = Parameters(θτRepository, pMax, vacancyRecoverDistance, typeDict)


# Initialize
simulator = Simulator(primaryVectors, boxSizes, inputGridVectors, latticeRanges, basis, basisTypes, parameters)
Save!(simulator)

function Irradiation(simulator::Simulator, energy::Float64)
    Restore!(simulator)
    simulator.nIrradiation += 1
    ionPosition = RandomInAnUnitGrapheneCell(1.39667) + [18.855045, 22.981482313368623, 20]
    ion = Atom(3, ionPosition, parameters)
    SetVelocityDirection!(ion, [0.,0.,-1.])
    SetEnergy!(ion,energy)
    push!(simulator, ion)
    Cascade!(ion, simulator)    
    BVacancy, NVacancy = CountVacancy(simulator)
    return BVacancy, NVacancy
end


Random.seed!(42)
computerNumberPerEnergy = 10000
vacancy_data = Dict{Int64, Vector{Tuple{Int64, Int64}}}()
x1 = [10*1.2^x for x in 0:35]
x2 = [x1[end]*1.7^x for x in 1:15]
energys = [x1..., x2...]
n = 0
for energy in energys
    global n += 1
    vacancyNumbers = Vector{Tuple{Int64, Int64}}(undef, computerNumberPerEnergy)
    @showprogress desc="In irradiation of energy order: $(n) ($(round(energy, digits=2)) eV)" for i in 1:computerNumberPerEnergy
        BVacancy, NVacancy = Irradiation(simulator, energy)
        vacancyNumbers[i] = (BVacancy, NVacancy)    
    end
    vacancy_data[n] = vacancyNumbers
end


log_info("Outputting data...")
open("vacancy.csv", "w") do f
    write(f, "n,energy,BVacancy,NVacancy\n")
    for n in sort(collect(keys(vacancy_data)))
        for (BVacancy, NVacancy) in vacancy_data[n]
            write(f, "$(n),$(energys[n]),$(BVacancy),$(NVacancy)\n")
        end
    end
end


