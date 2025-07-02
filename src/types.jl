#using PyCall

mutable struct Box
    vectors::Matrix{Float64}
    reciprocalVectors::Matrix{Float64}
    isOrthogonal::Bool
end


mutable struct Atom
    index::Int64  # never change
    isAlive::Bool
    type::Int64
    coordinate::Vector{Float64}
    cellIndex::Tuple{Int64, Int64, Int64}
    radius::Float64
    mass::Float64
    velocityDirection::Vector{Float64}
    energy::Float64
    Z::Float64

    dte::Float64
    bde::Float64

    numberOfEmptyCells::Int64

    # for atom_t
    pValue::Float64
    pPoint::Vector{Float64}
    pVector::Vector{Float64}
    pL::Float64

    # for atom_p
    lastTargets::Vector{Int64}


    latticePointIndex::Int64 # -1 for off lattice
    
    # for KMC 
    frequency::Float64 # Hz, s^-1
    frequencies::Vector{Float64} 
    finalLatticePointIndexs::Vector{Int64}
    eventIndex::Int64

    # for dynamic load 
    isNewlyLoaded::Bool
    latticeCoordinate::Vector{Float64}

end




mutable struct LatticePoint
    index::Int64
    type::Int64  # Initial Type, will not change
    coordinate::Vector{Float64}
    cellIndex::Tuple{Int64, Int64, Int64}
    environment::Vector{Int64}

    atomIndex::Int64 # -1 for vacancy
end


struct NeighborCellInfo
    index::NTuple{3, Int64}
    cross::NTuple{3, Int8} # 0 for no cross, 1 for hi, -1 for lo, eg. (0,0,1) for top 
end


mutable struct Cell
    # only for orthogonal box
    index::Tuple{Int64, Int64, Int64}
    atoms::Vector{Atom}
    latticePoints::Vector{LatticePoint}  
    ranges::Matrix{Float64}
    neighborCellsInfo::Array{NeighborCellInfo, 3}
    isExplored::Bool
    atomicDensity::Float64
    # for dynamic load
    latticeAtoms::Vector{Atom}
    isLoaded::Bool
    vacancies::Vector{Atom}
    isSavedLatticeRange::Bool
    latticeRanges::Matrix{Int64}
    isPushedNeighbor::Bool
end


function Cell(
    index::Tuple{Int64, Int64, Int64},
    atoms::Vector{Atom},
    latticePoints::Vector{LatticePoint},
    ranges::Matrix{Float64},
    #neighborCellsInfo::Dict{Vector{Int8}, NeighborCellInfo},
    neighborCellsInfo::Array{NeighborCellInfo, 3},
    isExplored::Bool,
    atomicDensity::Float64)           
    latticeAtoms = Vector{Atom}()
    isLoaded = false
    vacancies = Vector{Atom}()
    isSavedLatticeRange = false
    latticeRanges = Matrix{Int64}(undef, 3, 2)
    isPushedNeighbor = false
    return Cell(index, atoms, latticePoints, ranges, neighborCellsInfo, isExplored, atomicDensity, 
                    latticeAtoms, isLoaded, vacancies, isSavedLatticeRange, latticeRanges, isPushedNeighbor)     
end


macro cell_storage_type()
    if IS_DYNAMIC_LOAD
        return :(Dict{Tuple{Int64, Int64, Int64}, Cell}) #:(SparseVector{Cell})
    else
        return :(Array{Cell, 3})
    end
end

mutable struct Grid
    cells::@cell_storage_type()
    vectors::Matrix{Float64}
    sizes::Vector{Int64}      
    cellVolume::Float64
end 


struct ConstantsByType
    V_upterm::Dict{Vector{Int64}, Float64}
    a_U::Dict{Vector{Int64}, Float64}
    E_m::Dict{Int64, Float64}
    S_e_upTerm::Dict{Vector{Int64}, Float64}
    S_e_downTerm::Dict{Vector{Int64}, Float64}
    x_nl::Dict{Vector{Int64}, Float64}
    a::Dict{Vector{Int64}, Float64}
    Q_nl::Dict{Vector{Int64}, Float64}  
    Q_loc::Dict{Vector{Int64}, Float64}
    qMax::Dict{Vector{Int64}, Float64}
    sigma::Dict{Int64, Float64}
end


struct Element
    name::String
    radius::Float64
    mass::Float64
    Z::Float64
    dte::Float64
    bde::Float64
    alpha::Float64
    beta::Float64
end


mutable struct Parameters
    primaryVectors::Matrix{Float64}
    primaryVectors_INV::Matrix{Float64} # not a input 
    latticeRanges::Matrix{Int64}
    basisTypes::Vector{Int64}
    basis::Matrix{Float64}
    θτRepository::String
    pMax::Float64
    vacancyRecoverDistance_squared::Float64
    typeDict::Dict{Int64, Element}
    #optional 
    periodic::Vector{Bool}
    isOrthogonal::Bool
    isPrimaryVectorOrthogonal::Bool  # not a input 
    EPowerRange::StepRangeLen{Float64, Base.TwicePrecision{Float64}, Base.TwicePrecision{Float64}, Int64}
    pRange::StepRangeLen{Float64, Base.TwicePrecision{Float64}, Base.TwicePrecision{Float64}, Int64}
    stopEnergy::Float64
    DebyeTemperature::Float64
    pLMax::Float64     
    isDumpInCascade::Bool
    dumpFolder::String
    isLog::Bool
    DTEMode::Int64 
    #soapParameters::Vector{Float64}
    DTEFile::String
    isKMC::Bool
    nu_0_dict::Dict{Int64, Float64}
    temperature::Float64
    temperature_kb::Float64
    perfectEnvIndex::Int64
    irrdiationFrequency::Float64
    nCascadeEveryLoad::Int64
    maxRSS::Int64
    isAmorphous::Bool
end


function Parameters(
    # required
    primaryVectors::Matrix{Float64},
    latticeRanges::Matrix{Int64},
    basisTypes::Vector{Int64},
    basis::Matrix{Float64},
    θτRepository::String,
    pMax::Float64,  
    vacancyRecoverDistance::Float64, 
    typeDict::Dict{Int64, Element};
    # optional
    periodic::Vector{Bool} = [true, true, false],
    isOrthogonal::Bool = true,
    EPowerRange::StepRangeLen{Float64, Base.TwicePrecision{Float64}, Base.TwicePrecision{Float64}, Int64} = -1.0:0.045:8.0,
    pRange::StepRangeLen{Float64, Base.TwicePrecision{Float64}, Base.TwicePrecision{Float64}, Int64} = 0.0:0.01:5.0,
    stopEnergy::Float64 = 10.0, 
    DebyeTemperature::Float64 = 1000.0, 
    pLMax::Float64 = 2.0, 
    isDumpInCascade::Bool = false, 
    dumpFolder::String = ".",
    isLog::Bool = false,
    DTEMode::Int64 = 1,
    #soapParameters::Vector{Float64} = [2.6, 8.0, 6.0],
    DTEFile::String="",
    isKMC::Bool = false,
    nu_0_dict::Dict{Int64, Float64} = Dict{Int64, Float64}(), # Hz, s^-1
    temperature::Float64 = 0.0,   # K
    perfectEnvIndex::Int64 = 0,
    irrdiationFrequency::Float64 = 0.0,
    nCascadeEveryLoad = 1,
    maxRSS::Int = 10,
    isAmorphous = false) # unit: GB
    temperature_kb = temperature * 8.61733362E-5 # eV
    primaryVectors_INV = inv(primaryVectors)
    if !isdir(θτRepository)
        error("θτRepository $(θτRepository) does not exist.")
    end
    isPrimaryVectorOrthogonal = (primaryVectors[1,2] == 0.0 && primaryVectors[1,3] == 0.0 && 
                    primaryVectors[2,1] == 0.0 && primaryVectors[2,3] == 0.0 && 
                    primaryVectors[3,1] == 0.0 && primaryVectors[3,2] == 0.0)
    vacancyRecoverDistance_squared = vacancyRecoverDistance * vacancyRecoverDistance
    maxRSS *= 1048576  # unit: kB
    return Parameters(primaryVectors, primaryVectors_INV, latticeRanges, basisTypes, basis,
                      θτRepository, pMax,  vacancyRecoverDistance_squared, typeDict,
                      periodic, isOrthogonal, isPrimaryVectorOrthogonal,
                      EPowerRange, pRange, stopEnergy, DebyeTemperature, pLMax, isDumpInCascade, dumpFolder, isLog,
                      DTEMode, 
                      #soapParameters, 
                      DTEFile,
                      isKMC, nu_0_dict, temperature, temperature_kb, perfectEnvIndex, irrdiationFrequency,
                      nCascadeEveryLoad, maxRSS, isAmorphous)
end 

mutable struct CollisionParamsBuffers
    tanφList::Vector{Float64}  
    tanψList::Vector{Float64}
    E_tList::Vector{Float64}
    x_pList::Vector{Float64}
    x_tList::Vector{Float64}
    Q_locList::Vector{Float64}
    function CollisionParamsBuffers()
        return new(
            Vector{Float64}(), Vector{Float64}(), Vector{Float64}(),
            Vector{Float64}(), Vector{Float64}(), Vector{Float64}()
        )
    end
end

mutable struct WorkBuffers
    coordinate::Vector{Float64}
    candidateTargets::Vector{Atom}
    collisionParames::CollisionParamsBuffers
    threadCandidates::Vector{Vector{Atom}}
    function WorkBuffers(max_threads::Int64=Threads.nthreads())
        candidateTargets = Vector{Atom}()
        sizehint!(candidateTargets, 100)
        collisionParams = CollisionParamsBuffers()
        threadCandidates = [Vector{Atom}() for _ in 1:max_threads]
        for tc in threadCandidates
            sizehint!(tc, 50)
        end
        return new(Vector{Float64}(undef, 3), candidateTargets, 
                  collisionParams, threadCandidates)
    end
end

function EnsureCollisionCapacity!(buffers::CollisionParamsBuffers, n::Int)
    if length(buffers.tanφList) < n
        resize!(buffers.tanφList, n)
        resize!(buffers.tanψList, n)
        resize!(buffers.E_tList, n)
        resize!(buffers.x_pList, n)
        resize!(buffers.x_tList, n)
        resize!(buffers.Q_locList, n)
    end
end

function ClearBuffers!(buffers::WorkBuffers)
    empty!(buffers.coordinate)
    empty!(buffers.targets)
    empty!(buffers.candidateTargets)
    for tc in buffers.threadCandidates
        empty!(tc)
    end
end

mutable struct Simulator
    atoms::Vector{Atom}
    latticePoints::Vector{LatticePoint}
    box::Box
    grid::Grid
    maxAtomID::Int64
    numberOfAtoms::Int64
    constantsByType::ConstantsByType
    isStore::Bool
    displacedAtoms::Vector{Int64}
    numberOfAtomsWhenStored::Int64
    nCascade::Int64
    nCollisionEvent::Int64
    exploredCells::Vector{Cell}
    θFunctions::Dict{Vector{Int64}, Function}
    τFunctions::Dict{Vector{Int64}, Function}
    #soap::PyObject
    environmentCut::Float64
    DTEData::Vector{Vector{Float64}}
    #for kmc 
    time::Float64
    frequency::Float64
    frequencies::Vector{Float64}
    mobileAtoms::Vector{Atom}
    #for dynamic load
    loadedCells::Vector{Cell}
    vacancies::Vector{Atom}
    numberOfVacancies::Int64
    maxVacancyID::Int64
    minLatticeAtomID::Int64
    # for debug
    debugAtom::Atom
    parameters::Parameters
    workBuffers::WorkBuffers
end



function Simulator(box::Box, inputGridVectors::Matrix{Float64}, parameters::Parameters)
    grid = CreateGrid(box, inputGridVectors)
    constantsByType = InitConstantsByType(parameters.typeDict, parameters)
    θFunctions, τFunctions = InitθτFunctions(parameters, constantsByType)
    #soap = InitSoap(parameters)
    if parameters.DTEMode == 2
        environmentCut, DTEData = LoadDTEData(parameters)
    else
        environmentCut, DTEData = -1.0, Vector{Vector{Float64}}()
    end
    time = 0.0
    frequency = 0.0
    frequencies = Vector{Float64}()
    mobileAtoms = Vector{Atom}()
    loadedCells = Vector{Cell}()
    vacancies = Vector{Atom}()
    nCollisionEvent = 0
    numberOfVacancies = 0
    maxVacancyID = 1E6
    minLatticeAtomID = 0
    debugAtom = Atom(1, [0.0,0.0,0.0], parameters)
    workBuffers = WorkBuffers()
    return Simulator(Vector{Atom}(), Vector{LatticePoint}(), 
                     box, grid, 
                     0, 0, 
                     constantsByType,
                     false, Vector{Int64}(), 0, 
                     0,nCollisionEvent,
                     Vector{Cell}(),
                     θFunctions, τFunctions,
                     #soap, 
                     environmentCut, DTEData, 
                     time, frequency, frequencies, mobileAtoms,
                     loadedCells, vacancies, numberOfVacancies, maxVacancyID,minLatticeAtomID,
                     debugAtom,
                     parameters,
                     workBuffers)  
end

