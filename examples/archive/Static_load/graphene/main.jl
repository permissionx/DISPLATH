const IS_DYNAMIC_LOAD = false
home = ENV["ARCS_HOME"]
include(home * "/src/DISPLATH.jl")
seed = 42; const THREAD_RNG = [StableRNG(seed + t) for t in 1:Threads.nthreads()]


# Parameters
pMax = 1.45
vacancyRecoverDistance = 1.3
parameters = Parameters(pMax, vacancyRecoverDistance;
                        stopEnergy=0.1)


# Material
a = 1.42
b = 3.35
primaryVectors = [3.0*a 0.0 0.0; 0.0 3.0^0.5*a 0.0; 0.0 0.0 b]
latticeRanges = [0 10; 0 20; 5 6]   
basis = [0.0 0.0 0.0; 1.0/3.0 0.0 0.0; 1.0/2.0 1.0/2.0 0.0; 5.0/6.0 1.0/2.0 0.0]
basisTypes = [1, 1, 1, 1]
typeDict = Dict(
    1 => Element("C", 22.0, 11.0),   # dte, binding energy 
    2 => Element("Ne", 0.1, 0.1)
)
boxSizes = [10, 20, 10]
inputGridVectors = [2.1 0.0 0.0; 0.0 a*2.1 0.0; 0.0 0.0 a*2.1]  # never be same as primaryVectors 
material = Material(primaryVectors, latticeRanges, basisTypes, basis, typeDict,
                    boxSizes, inputGridVectors, parameters)
# Process 
simulator = Simulator(material, parameters)  
Save!(simulator)  
@dump "init.dump" simulator.atoms 




@showprogress for energy in 100.0:100.0:1000.0
    mean_nV = 0.0
    for i in 1:10000
        Restore!(simulator)
        ionPosition =  RandomInSquare(41.6, 48.19) + [0.1, 0.1, 33.0]
        Irradiation!(simulator, energy, ionPosition, [0.0,0.0,-1.0], 2, parameters)
        _, Vs = DefectStatics(simulator)
        nV = length(Vs)
        mean_nV += nV
    end
    mean_nV /= 10000
    @dump "final.dump" simulator.atoms 
    @record "nV.csv" "$(energy),$(mean_nV)" "energy,nV"
end

