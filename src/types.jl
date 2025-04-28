using PyCall

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
    cellIndex::Vector{Int64}
    radius::Float64
    mass::Float64
    velocityDirection::Vector{Float64}
    energy::Float64
    Z::Float64

    dte::Float64
    bde::Float64

    # for atom_t
    pValue::Dict{Int64, Float64}
    pPoint::Dict{Int64, Vector{Float64}}
    pVector::Dict{Int64, Vector{Float64}}
    pL::Dict{Int64, Float64}

    # for atom_p
    lastTargets::Vector{Int64}


    latticePointIndex::Int64 # -1 for off lattice
    
    # for KMC 
    probability::Float64
    probabilities::Vector{Float64} 
    finalLatticePointIndexs::Vector{Int64}
end


mutable struct LatticePoint
    index::Int64
    type::Int64  # Initial Type, will not change
    coordinate::Vector{Float64}
    cellIndex::Vector{Int64}
    environment::Vector{Int64}

    atomIndex::Int64 # -1 for vacancy
end


struct NeighborCellInfo
    index::Vector{Int64}
    cross::Vector{Int64} # 0 for no cross, 1 for hi, -1 for lo, eg. [0,0,1] for top 
end



mutable struct GridCell
    # only for orthogonal box
    index::Vector{Int64}
    atoms::Vector{Int64}  
    latticePoints::Vector{Int64}  
    ranges::Matrix{Float64}
    centerCoordinate::Vector{Float64}
    neighborCellsInfo::Dict{Vector{Int64}, NeighborCellInfo}
    isExplored::Bool
    atomicDensity::Float64
end


mutable struct CellGrid
    cells::Array{GridCell, 3}
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
    qMax_squared::Dict{Vector{Int64}, Float64}
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


struct Parameters
    θτRepository::String
    pMax::Float64
    vacancyRecoverDistance_squared::Float64
    typeDict::Dict{Int64, Element}
    #optional 
    periodic::Vector{Bool}
    isOrthogonal::Bool
    E_p_power_range::StepRangeLen{Float64, Base.TwicePrecision{Float64}, Base.TwicePrecision{Float64}, Int64}
    p_range::StepRangeLen{Float64, Base.TwicePrecision{Float64}, Base.TwicePrecision{Float64}, Int64}
    stopEnergy::Float64
    pLMax::Float64     
    isDumpInCascade::Bool
    isLog::Bool
    DTEMode::Int64 
    soapParameters::Vector{Float64}
    DTEFile::String
end


function Parameters(
    # required
    θτRepository::String,
    pMax::Float64,  
    vacancyRecoverDistance::Float64, 
    typeDict::Dict{Int64, Element};
    # optional
    periodic::Vector{Bool} = [true, true, false],
    isOrthogonal::Bool = true,
    E_p_power_range::StepRangeLen{Float64, Base.TwicePrecision{Float64}, Base.TwicePrecision{Float64}, Int64} = -1.0:0.045:8.0,
    p_range::StepRangeLen{Float64, Base.TwicePrecision{Float64}, Base.TwicePrecision{Float64}, Int64} = 0.0:0.01:2.0,
    stopEnergy::Float64 = 10.0, 
    pLMax::Float64 = 2.0, 
    isDumpInCascade::Bool = false, 
    isLog::Bool = false,
    DTEMode::Int64 = 1,
    soapParameters::Vector{Float64} = [2.6, 8.0, 6.0],
    DTEFile::String="")
    if !isdir(θτRepository)
        error("θτRepository $(θτRepository) does not exist.")
    end
    vacancyRecoverDistance_squared = vacancyRecoverDistance * vacancyRecoverDistance
    return Parameters(θτRepository, pMax,  vacancyRecoverDistance_squared, typeDict, 
                      periodic, isOrthogonal,
                      E_p_power_range, p_range, stopEnergy, pLMax, isDumpInCascade, isLog,
                      DTEMode, soapParameters, DTEFile)
end 


mutable struct Simulator
    atoms::Vector{Atom}
    latticePoints::Vector{LatticePoint}
    box::Box
    cellGrid::CellGrid
    maxAtomID::Int64
    numberOfAtoms::Int64
    constantsByType::ConstantsByType
    isStore::Bool
    displacedAtoms::Vector{Int64}
    atomNumberWhenStore::Int64
    nIrradiation::Int64
    exploredCells::Vector{GridCell}
    θFunctions::Dict{Vector{Int64}, Function}
    τFunctions::Dict{Vector{Int64}, Function}
    soap::PyObject
    environmentCut::Float64
    DTEData::Vector{Vector{Float64}}
    parameters::Parameters
    #for kmc 
    probability::Float64
    probabilities::Vector{Float64}
    mobileAtoms::Vector{Atom}
end

