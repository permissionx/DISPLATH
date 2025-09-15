home = ENV["ARCS_HOME"]
const BAlpha = 1.5
const IS_DYNAMIC_LOAD = false
include(home * "/src/DISPLATH.jl")

primaryVectors = [1.0 0. 0.; 0. 1. 0.; 0. 0. 1.]
latticeRanges = [0 10; 0 10; 0 10]
basisTypes = [1]
basis = [0.0 0.0 0.0]
pMax = 4.9
vacancyRecoverDistance = 0.0 
typeDict = Dict(
    1 => Element("C", 22.0, 11.0),  
    2 => Element("Ne", 0.1, 0.1)  
)
stopEnergy = 0.1 


parameters = Parameters(primaryVectors, latticeRanges, basisTypes, basis, pMax, vacancyRecoverDistance, typeDict;
                        stopEnergy=stopEnergy)

boxSizes = [10, 20, 10]
box = CreateBoxByPrimaryVectors(parameters.primaryVectors, boxSizes)
inputGridVectors = [3.4 0.0 0.0; 0.0 3.4 0.0; 0.0 0.0 3.4]
simulator = Simulator(box, inputGridVectors, parameters)


E = 500.
p = 1.716330802791826
atom_p = Atom(2, [5., 5., 6.], parameters)
SetVelocityDirection!(atom_p, [0.,0.,-1.])
SetEnergy!(atom_p, E)
push!(simulator, atom_p)

atom_t =  Atom(1, [5., 5.0-p , 4.], parameters)
push!(simulator, atom_t)

Cascade!(atom_p, simulator)


E_p = E
mass_p = 20.
mass_t = 12.
type_p = 2
type_t = 1
@show BCA.θτ(E_p, mass_p, mass_t, type_p, type_t, p, simulator.constantsByType)
