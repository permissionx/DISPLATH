home = ENV["ARCS_HOME"]
const IS_DYNAMIC_LOAD = false
include(home * "/src/DISPLATH.jl")

# Materials &  box
a = 1.42
b = 3.35
primaryVectors = [3.0*a 0.0 0.0; 0.0 3.0^0.5*a 0.0; 0.0 0.0 b]

boxSizes = [10, 20, 10]
periodic = [true, true, false]
latticeRanges = [0 10; 0 20; 5 6]   
basis = [0.0 0.0 0.0; 1.0/3.0 0.0 0.0; 1.0/2.0 1.0/2.0 0.0; 5.0/6.0 1.0/2.0 0.0]
basisTypes = [1, 1, 1, 1]

inputGridVectors = [2.1 0.0 0.0; 0.0 a*2.1 0.0; 0.0 0.0 a*2.1]  # never be same as primaryVectors 


typeDict = Dict(
    1 => Element("C", 22.0, 11.0),   # dte, binding energy 
    2 => Element("Ar", 0.1, 0.1)
)

# Parameters
pMax = 1.45  # From I2DM cut_radius approximation
vacancyRecoverDistance = 1.3  # From I2DM capture_radius setting
seed = 42; const THREAD_RNG = [StableRNG(seed + t) for t in 1:Threads.nthreads()]
stopEnergy = 0.1
parameters = Parameters(primaryVectors, latticeRanges, basisTypes, basis, pMax, vacancyRecoverDistance, typeDict;
                        stopEnergy=stopEnergy,periodic=periodic)


# Process
simulator = Simulator(boxSizes, inputGridVectors, parameters)  
@dump "init.dump" simulator.atoms 


function Irradiation(simulator::Simulator, energy::Float64)
    Restore!(simulator)
    ionPosition =  RandomInSquare(41.6, 48.19) + [0.1, 0.1, 33.0]
    ion = Atom(2, ionPosition, parameters)
    SetVelocityDirection!(ion, [0.0, 0.0, -1.0])
    SetEnergy!(ion,energy)
    push!(simulator, ion)
    Cascade!(ion, simulator)
    _, Vs = DefectStatics(simulator)
    return length(Vs)
end

en = 1:1:7
Save!(simulator)
@showprogress for ene in en
    energy = 10.0^ene
    mean_nV = 0.0
    mean_single = 0.0
    mean_double = 0.0
    for i in 1:10000
        nV = Irradiation(simulator::Simulator, energy::Float64)
        mean_nV += nV
        if nV <= 1 && nV > 0
            mean_single += 1
        elseif nV > 1
            mean_double += 1
        end
    end
    mean_nV /= 10000
    mean_single /= 10000
    mean_double /= 10000
    @dump "final.dump" simulator.atoms 
    @record "nV.csv" "$(energy),$(mean_nV),$(mean_single),$(mean_double)" "energy,nV,single,double"
end
