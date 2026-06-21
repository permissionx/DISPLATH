const IS_DYNAMIC_LOAD = false
home = ENV["ARCS_HOME"]
include(home * "/src/DISPLATH.jl")
seed = 42; const THREAD_RNG = [StableRNG(seed + t) for t in 1:Threads.nthreads()]


# Parameters
pMax = 1.45
vacancyRecoverDistance = 1.3
parameters = Parameters(pMax, vacancyRecoverDistance; isDumpInCascade=false,
                        stopEnergy=0.1)


# Material
a = 1.45
d = 3.34
Nlayers = 3   #  modify 
primaryVectors = [a*3 0.0 0.0; 0.0 3.0^0.5*a 0.0; 0.0 0.0 d*2]
latticeRanges = [0 10; 0 20; 2 2+Nlayers]     
#basis = [0.0 0.0 0.0; 1.0/3.0 0.0 0.0; 1.0/2.0 1.0/2.0 0.0; 5.0/6.0 1.0/2.0 0.0; 
#         1.0/6.0 1.0/2.0 0.5; 1.0/3.0+1.0/6.0 0.0+1.0/2.0 0.5; 1.0/2.0+1.0/6.0 1.0/2.0+1.0/2.0 0.5; 5.0/6.0+1.0/6.0 1.0/2.0+1.0/2.0 0.5;]
basis = [0.0 0.0 0.0; 1.0/3.0 0.0 0.0; 1.0/2.0 1.0/2.0 0.0; 5.0/6.0 1.0/2.0 0.0; 
         0.0 0.0 0.5; 1.0/3.0 0.0 0.5; 1.0/2.0 1.0/2.0 0.5; 5.0/6.0 1.0/2.0 0.5]
basisTypes = [1, 2, 1, 2, 2, 1, 2, 1]
typeDict = Dict(
    1 => Element("N", 7.9, 5.0),   # dte, binding energy 
    2 => Element("B", 10.1, 5.1), 
    3 => Element("O", 5.1, 5.1)
)
boxSizes = [10, 20, 10]
inputGridVectors = [2.1 0.0 0.0; 0.0 a*2.1 0.0; 0.0 0.0 a*2.1]  # never be same as primaryVectors 
material = Material(primaryVectors, latticeRanges, basisTypes, basis, typeDict,
                    boxSizes, inputGridVectors, parameters)
# Process 
simulator = Simulator(material, parameters)  
Save!(simulator)  
@dump "init.dump" simulator.atoms 


@showprogress for energy in 100.0:100.0:400.0
    Restore!(simulator)
    for i in 1:10
        ionPosition =  RandomInSquare(43.2, 49.0) + [0.1, 0.1, 43.0]
        Irradiation!(simulator, energy, ionPosition, [0.0,0.0,-1.0], 3, parameters)
    end
    _, Vs = DefectStatics(simulator)
    nV = length(Vs)
    @dump "final.dump" simulator.atoms 
    @record "nV.csv" "$(energy),$(nV)" "energy,nV"
end

