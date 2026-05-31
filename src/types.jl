#using PyCall

mutable struct Box
    vectors::Matrix{Float64}
    reciprocalVectors::Matrix{Float64}
    isOrthogonal::Bool
end


macro load_variant(dynamic_expr, static_expr)
    if !isdefined(__module__, :IS_DYNAMIC_LOAD)
        error("IS_DYNAMIC_LOAD must be defined before including DISPLATH types")
    end
    return getfield(__module__, :IS_DYNAMIC_LOAD) ? esc(dynamic_expr) : esc(static_expr)
end


@load_variant(
    begin
        mutable struct Atom
            index::Int64  # never change
            isAlive::Bool
            type::Int64
            coordinate::Vector{Float64}
            cellIndex::Tuple{Int64, Int64, Int64}
        end
    end,
    begin
    mutable struct Atom
        index::Int64  # never change
        isAlive::Bool
        type::Int64
        coordinate::Vector{Float64}
        cellIndex::Tuple{Int64, Int64, Int64}
        radius::Float64
        mass::Float64
        velocityDirection::SVector{3,Float64}
        energy::Float64
        Z::Float64

        dte::Float64
        bde::Float64

        emptyPath::Float64

        # for atom_t
        pValue::Float64
        pPoint::SVector{3,Float64}
        pVector::SVector{3,Float64}
        pL::Float64
        pAtomIndex::Int64
        pDirection::Vector{Float64}

        # for atom_p
        lastTargets::Vector{Int64}

        latticePointIndex::Int64 # -1 for off lattice

        # for KMC
        frequency::Float64
        frequencies::Vector{Float64}
        finalLatticePointIndexs::Vector{Int64}
        eventIndex::Int64

        # for dynamic load, kept here so static mode retains the main-branch shape
        isNewlyLoaded::Bool
        latticeCoordinate::SVector{3,Float64}
        indexInCell::Int64
    end
    end
)

struct Material
    box::Box
    atoms::Vector{Atom}
    inputGridVectors::Matrix{Float64}
end

mutable struct LatticePoint
    index::Int64
    type::Int64  # Initial Type, will not change
    coordinate::Vector{Float64}
    cellIndex::Tuple{Int64, Int64, Int64}
    environment::Vector{Int64}

    atomIndex::Int64 # -1 for vacancy
end


mutable struct NeighborCellInfo
    index::NTuple{3, Int64}
    cross::NTuple{3, Int8} # 0 for no cross, 1 for hi, -1 for lo, eg. (0,0,1) for top 
end

struct TargetCandidate
    index::Int64
    type::Int64
    cellIndex::Tuple{Int64, Int64, Int64}
    isLatticeAtom::Bool
    indexInCell::Int64
    coordinate::SVector{3,Float64}
    pValue::Float64
    pPoint::SVector{3,Float64}
    pVector::SVector{3,Float64}
    pL::Float64
end

struct AtomDynamics
    velocityDirection::SVector{3,Float64}
    energy::Float64
end

const ZERO_VELOCITY_DIRECTION = SVector{3,Float64}(0.0, 0.0, 0.0)
const ZERO_ATOM_DYNAMICS = AtomDynamics(ZERO_VELOCITY_DIRECTION, 0.0)


@load_variant(
    begin
        mutable struct Cell
            index::Tuple{Int64, Int64, Int64}
            atoms::Vector{Atom}
            vacancies::Vector{Atom}
        end
    end
,
    begin
    mutable struct Cell
        # only for orthogonal box
        index::Tuple{Int64, Int64, Int64}
        atoms::Vector{Atom}
        latticePoints::Vector{LatticePoint}
        ranges::Matrix{Float64}
        neighborCellsInfo::Array{NeighborCellInfo, 3}
        isExplored::Bool
        atomicDensity::Float64
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
    end
)

struct CellStd
    atoms::Vector{Atom}
    function CellStd()
        return new(Vector{Atom}())
    end
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
    V_upterm::Dict{Tuple{Int64, Int64}, Float64}
    a_U::Dict{Tuple{Int64, Int64}, Float64}
    E_m::Dict{Int64, Float64}
    S_e_upTerm::Dict{Tuple{Int64, Int64}, Float64}
    S_e_downTerm::Dict{Tuple{Int64, Int64}, Float64}
    x_nl::Dict{Tuple{Int64, Int64}, Float64}
    a::Dict{Tuple{Int64, Int64}, Float64}
    Q_nl::Dict{Tuple{Int64, Int64}, Float64}
    Q_loc::Dict{Tuple{Int64, Int64}, Float64}
    qMax::Dict{Tuple{Int64, Int64}, Float64}
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
    pMax_squared::Float64 # automatic
    vacancyRecoverDistance_squared::Float64
    typeDict::Dict{Int64, Element}
    #optional 
    periodic::Vector{Bool}
    isOrthogonal::Bool
    isPrimaryVectorOrthogonal::Bool  # not a input 
    EPowerRange::StepRangeLen{Float64, Base.TwicePrecision{Float64}, Base.TwicePrecision{Float64}, Int64}
    pPowerRange::StepRangeLen{Float64, Base.TwicePrecision{Float64}, Base.TwicePrecision{Float64}, Int64}
    stopEnergy::Float64
    isNonQnl::Bool
    DebyeTemperature::Float64
    isDumpInCascade::Bool
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
    amorphousLength::Float64
    amorphousHeight::Float64
    infiniteLength::Float64
    debugMode::Bool
end


function Parameters(
    # required
    primaryVectors::Matrix{Float64},
    latticeRanges::Matrix{Int64},
    basisTypes::Vector{Int64},
    basis::Matrix{Float64},
    pMax::Float64,  
    vacancyRecoverDistance::Float64, 
    typeDict::Dict{Int64, Element};
    # optional
    periodic::Vector{Bool} = [true, true, false],
    isOrthogonal::Bool = true,
    θτRepository::String = ENV["ARCS_REPO"] * "/thetatau_repository/",
    EPowerRange::StepRangeLen{Float64, Base.TwicePrecision{Float64}, Base.TwicePrecision{Float64}, Int64} = -1.0:0.045:8.0,
    pPowerRange::StepRangeLen{Float64, Base.TwicePrecision{Float64}, Base.TwicePrecision{Float64}, Int64} = -10.0:0.01:1.0,
    stopEnergy::Float64 = 0.1, 
    isNonQnl::Bool = false,     #only works when static loading  
    DebyeTemperature::Float64 = 519.0,  # K
    isDumpInCascade::Bool = false, 
    DTEMode::Int64 = 1,
    #soapParameters::Vector{Float64} = [2.6, 8.0, 6.0],
    DTEFile::String="",
    isKMC::Bool = false,
    nu_0_dict::Dict{Int64, Float64} = Dict{Int64, Float64}(), # Hz, s^-1
    temperature::Float64 = 0.0,   # K
    perfectEnvIndex::Int64 = 0,
    irrdiationFrequency::Float64 = 0.0,
    nCascadeEveryLoad = 100,
    maxRSS::Int = 20, # unit: GB
    isAmorphous::Bool = false,
    amorphousLength::Float64 = -100.0,
    infiniteLength::Float64 = 1000.0,
    debugMode::Bool = false) 

    pMax_squared = pMax * pMax 
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
    amorphousHeight =  latticeRanges[3,2] * primaryVectors[3,3] - amorphousLength
    return Parameters(primaryVectors, primaryVectors_INV, latticeRanges, basisTypes, basis,
                      θτRepository, pMax, pMax_squared, vacancyRecoverDistance_squared, typeDict,
                      periodic, isOrthogonal, isPrimaryVectorOrthogonal,
                      EPowerRange, pPowerRange, stopEnergy, isNonQnl, DebyeTemperature, isDumpInCascade, 
                      DTEMode, 
                      #soapParameters, 
                      DTEFile,
                      isKMC, nu_0_dict, temperature, temperature_kb, perfectEnvIndex, irrdiationFrequency,
                      nCascadeEveryLoad, maxRSS, isAmorphous, amorphousLength, amorphousHeight, infiniteLength,
                      debugMode)
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
    coordinates::Vector{Vector{Float64}}
    candidateTargets::Vector{TargetCandidate}
    collisionParames::CollisionParamsBuffers
    threadCandidates::Vector{Vector{TargetCandidate}}
    neighborCellsInfos::Vector{Array{NeighborCellInfo, 3}}
    latticeSiteCoordinates::Vector{Dict{Tuple{Int64, Int64, Int64}, Vector{SVector{3,Float64}}}}
    lastTargets::Dict{Int64, Vector{Int64}}
    atomDynamics::Dict{Int64, AtomDynamics}
    function WorkBuffers(max_threads::Int64=Threads.nthreads())
        coordinates = [Vector{Float64}(undef, 3) for _ in 1:max_threads]
        candidateTargets = Vector{TargetCandidate}()
        sizehint!(candidateTargets, 100)
        collisionParams = CollisionParamsBuffers()
        threadCandidates = [Vector{TargetCandidate}() for _ in 1:max_threads]
        for tc in threadCandidates
            sizehint!(tc, 50)
        end
        neighborCellsInfos = [_new_neighbor_info_buffer() for _ in 1:2]
        latticeSiteCoordinates = [Dict{Tuple{Int64, Int64, Int64}, Vector{SVector{3,Float64}}}() for _ in 1:max_threads]
        lastTargets = Dict{Int64, Vector{Int64}}()
        atomDynamics = Dict{Int64, AtomDynamics}()
        return new(coordinates, candidateTargets,
                  collisionParams, threadCandidates, neighborCellsInfos, latticeSiteCoordinates,
                  lastTargets, atomDynamics)
    end
end

function _new_neighbor_info_buffer()
    buffer = Array{NeighborCellInfo, 3}(undef, 3, 3, 3)
    for i in eachindex(buffer)
        buffer[i] = NeighborCellInfo((0, 0, 0), (Int8(0), Int8(0), Int8(0)))
    end
    return buffer
end

function LastTargets!(atom::Atom, simulator)
    return get!(simulator.workBuffers.lastTargets, atom.index) do
        Vector{Int64}()
    end
end

function ClearLastTargets!(atom::Atom, simulator)
    targets = get(simulator.workBuffers.lastTargets, atom.index, nothing)
    targets === nothing || empty!(targets)
    return nothing
end

function AtomDynamics!(atom::Atom, simulator)
    return AtomDynamics!(atom.index, simulator)
end

function AtomDynamics!(index::Int64, simulator)
    return get(simulator.workBuffers.atomDynamics, index, ZERO_ATOM_DYNAMICS)
end

function ClearAtomDynamics!(atom::Atom, simulator)
    return ClearAtomDynamics!(atom.index, simulator)
end

function ClearAtomDynamics!(index::Int64, simulator)
    delete!(simulator.workBuffers.atomDynamics, index)
    return nothing
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
    empty!(buffers.candidateTargets)
    for tc in buffers.threadCandidates
        empty!(tc)
    end
    ClearLatticeSiteCoordinateCaches!(buffers)
    empty!(buffers.lastTargets)
    empty!(buffers.atomDynamics)
end

function ClearLatticeSiteCoordinateCaches!(buffers::WorkBuffers)
    for cache in buffers.latticeSiteCoordinates
        empty!(cache)
    end
    return nothing
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
    θFunctions::Dict{Tuple{Int64, Int64}, Function}
    τFunctions::Dict{Tuple{Int64, Int64}, Function}
    uniformDensity::Float64
    #soap::PyObject
    environmentCut::Float64
    DTEData::Vector{Vector{Float64}}
    #for kmc 
    time::Float64
    frequency::Float64
    frequencies::Vector{Float64}
    mobileAtoms::Vector{Atom}
    #for dynamic load
    vacancies::Vector{Atom}
    numberOfVacancies::Int64
    maxVacancyID::Int64
    minLatticeAtomID::Int64
    deprecatedCellKeys::Set{Tuple{Int64, Int64, Int64}}
    preservedCellKeys::Set{Tuple{Int64, Int64, Int64}}
    attempedDeCellKeys::Set{Tuple{Int64, Int64, Int64}}
    cellStd::CellStd
    cellLatticeAtomNumber::Int64
    # for debug
    debugAtoms::Vector{Atom}
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
    vacancies = Vector{Atom}()
    nCollisionEvent = 0
    numberOfVacancies = 0
    maxVacancyID = 1E6
    minLatticeAtomID = 0
    debugAtoms = Atom[]
    workBuffers = WorkBuffers()
    uniformDensity = length(parameters.basisTypes) / (parameters.primaryVectors[1,1] * parameters.primaryVectors[2,2] * parameters.primaryVectors[3,3])
    cellStd = CellStd()
    deprecatedCellKeys = Set{Tuple{Int64, Int64, Int64}}()
    preservedCellKeys = Set{Tuple{Int64, Int64, Int64}}()
    attempedDeCellKeys = Set{Tuple{Int64, Int64, Int64}}()
    cellLatticeAtomNumber = 0
    return Simulator(Vector{Atom}(), Vector{LatticePoint}(), 
                     box, grid, 
                     0, 0, 
                     constantsByType,
                     false, Vector{Int64}(), 0, 
                     0,nCollisionEvent,
                     Vector{Cell}(),
                     θFunctions, τFunctions,
                     uniformDensity,
                     #soap, 
                     environmentCut, DTEData, 
                     time, frequency, frequencies, mobileAtoms,
                     vacancies, numberOfVacancies, maxVacancyID, minLatticeAtomID, 
                     deprecatedCellKeys, preservedCellKeys, attempedDeCellKeys, 
                     cellStd, cellLatticeAtomNumber, 
                     debugAtoms,
                     parameters,
                     workBuffers)  
end
