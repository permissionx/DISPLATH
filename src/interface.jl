function Parameters(pMax::Float64, vacancyRecoverDistance::Float64; kwargs...)
    # non lattice info
    primaryVectors = [1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0]
    latticeRanges = [0 1; 0 1; 0 1]
    basis = [0.0 0.0 0.0]
    basisTypes = [1]  
    typeDict = Dict{Int64, Element}()
    parameters = Parameters(primaryVectors, latticeRanges, basisTypes, basis, pMax, vacancyRecoverDistance, typeDict; kwargs...)
    return parameters
end

function Material(
    primaryVectors::Matrix{Float64},
    latticeRanges::Matrix{Int64},
    basisTypes::Vector{Int64},
    basis::Matrix{Float64},
    typeDict::Dict{Int64, Element},
    boxSizes::Vector{Int64}, 
    inputGridVectors::Matrix{Float64},
    parameters::Parameters)
    atoms = _PassingLatticeParamtersAndCreateAtoms(primaryVectors, latticeRanges, basisTypes, basis, typeDict, parameters)
    box = CreateBoxByPrimaryVectors(parameters.primaryVectors, boxSizes)
    material = Material(box, atoms, inputGridVectors)
    return material
end

function Material(
    primaryVectors::Matrix{Float64},
    latticeRanges::Matrix{Int64},
    basisTypes::Vector{Int64},
    basis::Matrix{Float64},
    typeDict::Dict{Int64, Element},
    boxVectors::Matrix{Float64}, 
    inputGridVectors::Matrix{Float64},
    parameters::Parameters)
    atoms = _PassingLatticeParamtersAndCreateAtoms(primaryVectors, latticeRanges, basisTypes, basis, typeDict, parameters)
    box = Box(boxVectors)
    material = Material(box, atoms, inputGridVectors)
    return material
end

function Material(fileName::String, typeDict::Dict{Int64, Element}, inputGridVectors::Matrix{Float64}, parameters::Parameters; replicate::Vector{Int64} = [1,1,1])
    # In this mode, primaryVectors, latticeRanges, basisTypes, basis in parameters are expired.
    if IS_DYNAMIC_LOAD
        error("Simulator from date file is not supported in dynamic load mode.")
    end 
    parameters.typeDict = typeDict
    box, atoms = LoadAtomsAndBoxFromDataFile(fileName; replicate=replicate)
    material = Material(box, atoms, inputGridVectors)
    return material
end

function Simulator(material::Material, parameters::Parameters)
    simulator = Simulator(material.box, material.atoms, material.inputGridVectors, parameters)
    return simulator
end


# utils
function _PassingLatticeParamtersAndCreateAtoms(
    primaryVectors::Matrix{Float64},
    latticeRanges::Matrix{Int64},
    basisTypes::Vector{Int64},
    basis::Matrix{Float64},
    typeDict::Dict{Int64, Element},
    parameters::Parameters)
    parameters.primaryVectors = primaryVectors
    parameters.latticeRanges = latticeRanges
    parameters.basisTypes = basisTypes
    parameters.basis = basis
    parameters.typeDict = typeDict
    if !IS_DYNAMIC_LOAD
        atoms = CreateAtomsByPrimaryVectors(parameters)
    else
        atoms = Atom[]
    end
    return atoms
end

