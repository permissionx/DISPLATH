import Base: push!
import Base: delete!
using .BCA.ConstantFunctions
using LinearAlgebra
using Interpolations   
using Dates
using ProgressMeter
using PyCall
using QuadGK
# @pyimport dscribe.descriptors as descriptors

function Box(Vectors::Matrix{Float64})
    # need to improve to detact if it is orithogonal. 
    println("Box with size of $(Vectors[1,1]) $(Vectors[2,2]) $(Vectors[3,3]) created.")
    return Box(Vectors, inv(Vectors'), true)
end 

function CreateBoxByPrimaryVectors(primaryVectors::Matrix{Float64}, sizes::Vector{Int64})
    vectors = primaryVectors .* sizes
    return Box(vectors)
end 


function Atom(type::Int64, coordinate::Vector{Float64}, parameters::Parameters)
    index = 0
    isAlive = true
    cellIndex = Vector{Int64}(undef, 3)
    velocityDirection = Float64[0.0,0.0,0.0]  # length： 0 or one
    energy = 0.0
    radius, mass, Z, dte, bde, _, _ = TypeToProperties(type, parameters.typeDict)
    pValue = Dict{Int64, Float64}() # key is the index of the atom_p
    pVector = Dict{Int64, Vector{Float64}}()
    pPoint = Dict{Int64, Vector{Float64}}()  
    lastTargets = Vector{Int64}()
    pL = Dict{Int64, Float64}()
    latticePointIndex = -1
    frequency = 0.0
    frequencies = Vector{Float64}()
    finalLatticePointEnvIndexs = Vector{Int64}()
    eventIndex = -1
    isNewlyLoaded = false
    lattcieCoordinate = Vector{Float64}(undef, 3)
    return Atom(index, isAlive, type, coordinate, cellIndex, 
                radius, mass, velocityDirection, energy, Z, 
                dte, bde,
                pValue, pVector, pPoint, pL, lastTargets,
                latticePointIndex,
                frequency, frequencies, finalLatticePointEnvIndexs, eventIndex, 
                isNewlyLoaded, lattcieCoordinate)
end


function TypeToProperties(type::Int64, typeDict::Dict{Int64, Element})
    if haskey(typeDict, type)
        element = typeDict[type]
        return element.radius, element.mass, element.Z, element.dte, element.bde, element.alpha, element.beta 
    else
        error("Unknown atom type: $type")
    end 
end 


function IterPushCellNeighbors!(cellGrid::CellGrid, gridCell::GridCell)
    # Direct triple loop implementation - much faster than recursion
    for delta_x in [-1, 0, 1]
        for delta_y in [-1, 0, 1]
            for delta_z in [-1, 0, 1]
                neighborKeys = Int8[delta_x, delta_y, delta_z]
                neighborIndex = Vector{Int64}(undef, 3)
                neighborCross = Vector{Int8}(undef, 3)
                
                # Calculate neighbor cell index and cross flags for each dimension
                for d in 1:3
                    delta = neighborKeys[d]
                    index = gridCell.index[d] + delta
                    cross = 0
                    if index < 1
                        index += cellGrid.sizes[d]
                        cross = -1
                    elseif index > cellGrid.sizes[d]
                        index -= cellGrid.sizes[d]
                        cross = 1
                    end
                    neighborIndex[d] = index
                    neighborCross[d] = cross
                end
                
                neighborCellInfo = NeighborCellInfo(neighborIndex, neighborCross)
                gridCell.neighborCellsInfo[neighborKeys] = neighborCellInfo
            end
        end
    end
end

function CreateCellGrid(box::Box, inputVectors::Matrix{Float64}, parameters::Parameters)
    if !box.isOrthogonal
        error("The box is not orthogonal, please use the orthogonal box.")
    end
    sizes = Vector{Int64}(undef, 3)
    vectors = Matrix{Float64}(undef, 3, 3)
    for d in 1:3
        sizes[d] = Int64(floor(box.vectors[d,d] / inputVectors[d,d]))
        vectors[d,d] = box.vectors[d,d] / sizes[d]
    end
    println("Cell grid: $(sizes[1]) x $(sizes[2]) x $(sizes[3]) = $(sizes[1]*sizes[2]*sizes[3]) cells with each cell size of $(vectors[1,1]) $(vectors[2,2]) $(vectors[3,3]).")
    
    cells = Array{GridCell, 3}(undef, sizes[1], sizes[2], sizes[3])
    @showprogress desc="Creating cell grid: " for x in 1:sizes[1]
        for y in 1:sizes[2]
            for z in 1:sizes[3]
                ranges = Matrix{Float64}(undef, 3, 2)
                ranges[1,1] = (x-1) * vectors[1,1]
                ranges[1,2] = x * vectors[1,1]
                ranges[2,1] = (y-1) * vectors[2,2]
                ranges[2,2] = y * vectors[2,2]
                ranges[3,1] = (z-1) * vectors[3,3]
                ranges[3,2] = z * vectors[3,3]  
                centerCoordinate = Vector{Float64}(undef, 3)  
                for d in 1:3
                    centerCoordinate[d] = (ranges[d,1] + ranges[d,2])/2
                end
                cells[x, y, z] = GridCell(Vector{Int64}([x,y,z]), Vector{Atom}(), Vector{LatticePoint}(), 
                                          ranges, centerCoordinate, 
                                          Dict{Vector{Int8}, NeighborCellInfo}(), false, 0.0)
            end
        end    
    end
    cellVolume = vectors[1,1] * vectors[2,2] * vectors[3,3]
    cellGrid = CellGrid(cells, vectors, sizes, cellVolume) 
    if ! parameters.isDynamicLoad
        @showprogress desc="Pushing cell neighbors: " for cell in cellGrid.cells
            IterPushCellNeighbors!(cellGrid, cell)
        end
    end
    println("Cell grid created.")
    return cellGrid
end





function InitConstantsByType(typeDict::Dict{Int64, Element}, parameters::Parameters) 
    V_upterm = Dict{Vector{Int64}, Float64}()
    a_U = Dict{Vector{Int64}, Float64}()
    E_m = Dict{Int64, Float64}()
    S_e_upTerm = Dict{Vector{Int64}, Float64}()
    S_e_downTerm = Dict{Vector{Int64}, Float64}()
    x_nl = Dict{Vector{Int64}, Float64}()
    a = Dict{Vector{Int64}, Float64}()
    Q_nl = Dict{Vector{Int64}, Float64}()  
    Q_loc = Dict{Vector{Int64}, Float64}()
    types = keys(typeDict)
    qMax = Dict{Vector{Int64}, Float64}()
    sigma = Dict{Int64, Float64}()
    for p in types
        for t in types
            radius_p, mass_p, Z_p, _, _, α_p, β_p = TypeToProperties(p, typeDict)
            radius_t, _, Z_t, _, _, _, _ = TypeToProperties(t, typeDict)
            V_upterm[[p,t]] = BCA.ConstantFunctions.V_upterm(Z_p, Z_t)
            a_U[[p,t]] = BCA.ConstantFunctions.a_U(Z_p, Z_t)
            E_m[p] = BCA.ConstantFunctions.E_m(Z_p, mass_p)
            S_e_upTerm[[p,t]] = BCA.ConstantFunctions.S_e_upTerm(p, Z_p, Z_t, mass_p, α_p)
            x_nl[[p,t]] = BCA.ConstantFunctions.x_nl(p, Z_p, Z_t, β_p)
            a[[p,t]] = BCA.ConstantFunctions.a(Z_p, Z_t)
            Q_nl[[p,t]] = BCA.ConstantFunctions.Q_nl(Z_p, Z_t, parameters.pMax)
            Q_loc[[p,t]] = BCA.ConstantFunctions.Q_loc(Z_p, Z_t)
            qMax[[p,t]] = radius_p + radius_t
            sigma[p] = TemperatureToSigma(parameters.temperature, parameters.DebyeTemperature, mass_p)
        end
    end
    return ConstantsByType(V_upterm, a_U, E_m, S_e_upTerm, S_e_downTerm, x_nl, a, Q_nl, Q_loc, qMax, sigma)
end


function InitθτFunctions(parameters::Parameters, constantsByType::ConstantsByType)
    typeDict = parameters.typeDict
    θFunctions = Dict{Vector{Int64}, Function}()
    τFunctions = Dict{Vector{Int64}, Function}()
    
    for type_p in keys(typeDict)
        for type_t in keys(typeDict)
            mass_p = typeDict[type_p].mass
        mass_t = typeDict[type_t].mass
        θInterpolation, τInterpolation = θτFunctions(mass_p, mass_t, type_p, type_t, constantsByType, parameters)
        
        θFunctions[[type_p, type_t]] = (E_p, p) -> θInterpolation(E_p, p)
        τFunctions[[type_p, type_t]] = (E_p, p) -> τInterpolation(E_p, p)
        end
    end
    println("All θ and τ functions initialized.\n")
    return θFunctions, τFunctions
end


function θτFunctions(mass_p::Float64, mass_t::Float64, type_p::Int64, type_t::Int64, constantsByType::ConstantsByType, parameters::Parameters)
    E_p_axis = Float64[]
    p_axis = Float64[]
    θMatrix = Matrix{Float64}(undef, 0, 0)
    τMatrix = Matrix{Float64}(undef, 0, 0)
    try
        E_p_axis, p_axis, θMatrix, τMatrix = LoadθτData(type_p, type_t, parameters)
        println("θ and τ functions for $(parameters.typeDict[type_p].name) to $(parameters.typeDict[type_t].name) Loaded.")
    catch
        println("Loading θ and τ data for $(parameters.typeDict[type_p].name) to $(parameters.typeDict[type_t].name) failed.")
        EPowerRange = parameters.EPowerRange    
        pRange = parameters.pRange
        nE = length(EPowerRange)
        np = length(pRange)
        θMatrix = Array{Float64, 2}(undef, nE, np)
        τMatrix = Array{Float64, 2}(undef, nE, np)
        @showprogress desc="Waiting for generating θ and τ data: " for (i, E_p_power) in enumerate(EPowerRange)
            E_p = 10.0^E_p_power
            for (j, p) in enumerate(pRange)
                θ, τ = BCA.θτ(E_p, mass_p, mass_t, type_p, type_t, p, constantsByType)
                θMatrix[i, j] = θ
                τMatrix[i, j] = τ
            end
        end
        E_p_axis = [10.0^E_p_power for E_p_power in EPowerRange]
        p_axis = collect(pRange)    
        SaveθτData(type_p, type_t, θMatrix, τMatrix, E_p_axis, p_axis, parameters)
    end

    # interpolate
    θFunction = interpolate((E_p_axis, p_axis), θMatrix, Gridded(Linear()))
    τFunction = interpolate((E_p_axis, p_axis), τMatrix, Gridded(Linear()))
    return θFunction, τFunction
end

function Simulator_dynamicLoad(boxSizes::Vector{Int64}, inputGridVectors::Matrix{Float64}, parameters::Parameters)
    println("Initializing the simulator...")
    box = CreateBoxByPrimaryVectors(parameters.primaryVectors, boxSizes)
    ranges = parameters.latticeRanges 
    atomNumber = (ranges[1,2] - ranges[1,1]) * (ranges[2,2] - ranges[2,1]) * (ranges[3,2] - ranges[3,1]) * length(parameters.basisTypes)
    println("$(atomNumber) atoms in the box.") 
    simulator = Simulator(box, inputGridVectors, parameters)
    vectors = parameters.primaryVectors     
    unitCellVolume = abs(dot(cross(vectors[:,1], vectors[:,2]), vectors[:,3]))
    for cell in simulator.cellGrid.cells
        cell.atomicDensity = length(parameters.basisTypes) / unitCellVolume
    end
    println("Simulator initialized.\n")
    return simulator    
end 

function Simulator(boxSizes::Vector{Int64}, inputGridVectors::Matrix{Float64}, parameters::Parameters)
    println("Initializing the simulator...")
    box = CreateBoxByPrimaryVectors(parameters.primaryVectors, boxSizes)
    simulator = Simulator(box, inputGridVectors, parameters)
    _initSimulatorAtoms!(simulator, parameters)
    return simulator    
end 

function Simulator(boxVectors::Matrix{Float64}, inputGridVectors::Matrix{Float64}, parameters::Parameters)
    println("Initializing the simulator...")
    box = Box(boxVectors)
    simulator = Simulator(box, inputGridVectors, parameters)
    _initSimulatorAtoms!(simulator, parameters)
    return simulator
end 

function _initSimulatorAtoms!(simulator::Simulator, parameters::Parameters)
    primaryVectors = parameters.primaryVectors
    latticeRanges = parameters.latticeRanges
    basis = parameters.basis
    basisTypes = parameters.basisTypes
    atomNumber = (latticeRanges[1,2] - latticeRanges[1,1]) * (latticeRanges[2,2] - latticeRanges[2,1]) * (latticeRanges[3,2] - latticeRanges[3,1]) * length(basisTypes)
    @showprogress desc="Creating atoms ($(atomNumber)): " for x in latticeRanges[1,1]:latticeRanges[1,2]-1
        for y in latticeRanges[2,1]:latticeRanges[2,2]-1    
            for z in latticeRanges[3,1]:latticeRanges[3,2]-1
                for i in 1:length(basisTypes)
                    reducedCoordinate = Float64[x,y,z] + basis[i, :]
                    coordinate = primaryVectors' * reducedCoordinate
                    atom = Atom(basisTypes[i], coordinate, parameters)
                    push!(simulator, atom)
                    latticePoint = LatticePoint(atom)
                    push!(simulator, latticePoint)
                end
            end
        end
    end
    println("$(simulator.numberOfAtoms) atoms created.")
    InitLatticePointEnvronment(simulator)
    for cell in simulator.cellGrid.cells
        cell.atomicDensity = length(cell.atoms) / simulator.cellGrid.cellVolume
    end 
    println("Simulator initialized.\n")
end

function LatticePoint(atom::Atom)
    environment = Vector{Int64}()
    return LatticePoint(copy(atom.index), copy(atom.type), 
                        copy(atom.coordinate), copy(atom.cellIndex), environment,
                        atom.index)
end

function WhichCell(coordinate::Vector{Float64}, cellGrid::CellGrid)
    cellIndex = Vector{Int64}(undef, 3)
    for d in 1:3
        cellIndex[d] = Int64(floor(coordinate[d] / cellGrid.vectors[d,d])) + 1
        if cellIndex[d] < 1 
            cellIndex[d] = 1
        elseif cellIndex[d] > cellGrid.sizes[d]
            cellIndex[d] = cellGrid.sizes[d]
        end
    end
    return cellIndex
end


function push!(simulator::Simulator, atom::Atom)
    atom.index = simulator.maxAtomID + 1
    simulator.maxAtomID += 1
    push!(simulator.atoms, atom)
    simulator.numberOfAtoms += 1
    cellIndex = WhichCell(atom.coordinate, simulator.cellGrid)
    atom.cellIndex .= cellIndex
    push!(simulator.cellGrid.cells[cellIndex[1], cellIndex[2], cellIndex[3]].atoms, atom)
end 


function push!(simulator::Simulator, latticePoint::LatticePoint)
    push!(simulator.latticePoints, latticePoint)
    push!(simulator.cellGrid.cells[latticePoint.cellIndex[1], latticePoint.cellIndex[2], latticePoint.cellIndex[3]].latticePoints, latticePoint)
    simulator.atoms[latticePoint.atomIndex].latticePointIndex = latticePoint.index
end 

function push!(cell::GridCell, atom::Atom, simulator::Simulator)
    atom.cellIndex .= cell.index
    push!(cell.atoms, atom)
    if !simulator.parameters.isDynamicLoad
        cell.atomicDensity = length(cell.atoms) / simulator.cellGrid.cellVolume
    end
end 


function delete!(simulator::Simulator, atom::Atom)
    originalCell = simulator.cellGrid.cells[atom.cellIndex[1], atom.cellIndex[2], atom.cellIndex[3]]
    deleteat!(originalCell.atoms, findfirst(a -> a.index == atom.index, originalCell.atoms))
    simulator.numberOfAtoms -= 1
    atom.isAlive = false 
    if atom.latticePointIndex != -1
        LeaveLatticePoint!(atom, simulator)
    end
end

function LeaveLatticePoint!(atom::Atom, simulator::Simulator; isUpdateEnv::Bool = true)
    AddToStore!(atom, simulator)       
    latticePoint = simulator.latticePoints[atom.latticePointIndex]
    latticePoint.atomIndex = -1
    atom.latticePointIndex = -1
    if isUpdateEnv && simulator.parameters.isKMC
        DeleteAtomEvents!(simulator, atom)
        UpdateEvents!(Set(latticePoint.environment), simulator)
    end
end


function delete!(cell::GridCell, atom::Atom, simulator::Simulator)
    if !atom.isAlive
        error("Atom $(atom.index) is not alive when deleting")
    end
    deleteat!(cell.atoms, findfirst(a -> a.index == atom.index, cell.atoms))
    atom.cellIndex = Vector{Int64}([-1, -1, -1])
    if !simulator.parameters.isDynamicLoad
        cell.atomicDensity = length(cell.atoms) / simulator.cellGrid.cellVolume  
    end
end


function DisplaceAtom!(atom::Atom, newPosition::Vector{Float64}, simulator::Simulator)
    if atom.type == 3 && newPosition[3] - atom.coordinate[3] < -10
        @show atom.coordinate
        @show newPosition
        @show simulator.nCollisionEvent
        println("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@")
    end
    if atom.type == 3
        Log("$(newPosition[3])\n",fileName="z.csv",type="a")
    end
    for d in 1:3
        # need to adappt non-periodic condition
        if newPosition[d] < 0
            if simulator.parameters.periodic[d] == false
                newPosition[d] = 0.01
            else
                newPosition[d] += simulator.box.vectors[d,d]
            end
        elseif newPosition[d] >= simulator.box.vectors[d,d]
            if simulator.parameters.periodic[d] == false
                newPosition[d] = simulator.box.vectors[d,d] - 0.01
            else
                newPosition[d] -= simulator.box.vectors[d,d]
            end
        end
    end
    SetCoordinate!(atom, newPosition)
    cellIndex = WhichCell(atom.coordinate, simulator.cellGrid)
    if cellIndex != atom.cellIndex
        ChangeCell!(atom, cellIndex, simulator)
    end
end


function ComputeDistance_squared(coordinate1::Vector{Float64}, coordinate2::Vector{Float64}, crossFlag::Vector{Int8}, box::Box)
    dv = VectorDifference(coordinate1, coordinate2, crossFlag, box)
    distance_squared = dv[1]* dv[1] + dv[2]*dv[2] + dv[3]  * dv[3]
    return distance_squared
end


function ComputeDistance(coordinate1::Vector{Float64}, coordinate2::Vector{Float64}, crossFlag::Vector{Int8}, box::Box)
    return sqrt(ComputeDistance_squared(coordinate1, coordinate2, crossFlag, box))
end


function ComputeVDistance(atom_p::Atom, atom_t::Atom, crossFlag::Vector{Int8}, box::Box)
    # v for atom_p
    dv = VectorDifference(atom_p.coordinate, atom_t.coordinate, crossFlag, box)
    return dv' * atom_p.velocityDirection
end


function VectorDifference(v1::Vector{Float64}, v2::Vector{Float64}, crossFlag::Vector{Int8}, box::Box)
    # return v2 - v1
    if crossFlag == Vector{Int64}([0,0,0])
        return v2 - v1
    end 
    result = Vector{Float64}(undef, 3)
    for d in 1:3
        dv = v2[d] - v1[d] + crossFlag[d] * box.vectors[d,d]
        result[d] = dv
    end
    return result
end


function ComputeP!(atom_p::Atom, atom_t::Atom, crossFlag::Vector{Int8}, box::Box)
    dv = VectorDifference(atom_p.coordinate, atom_t.coordinate, crossFlag, box)
    t = sum(dv .* atom_p.velocityDirection) / sum(atom_p.velocityDirection .* atom_p.velocityDirection)
    atom_t.pL[atom_p.index] = t
    atom_t.pPoint[atom_p.index] = atom_p.coordinate + t * atom_p.velocityDirection
    atom_t.pVector[atom_p.index] = atom_t.pPoint[atom_p.index] - atom_t.coordinate
    atom_t.pValue[atom_p.index] = norm(atom_t.pVector[atom_p.index])
end

function GetAllNeighbors(gridCell::GridCell, simulator::Simulator)
    cellGrid = simulator.cellGrid
    neighbors = Vector{Atom}()
    for (_, neighborCellInfo) in gridCell.neighborCellsInfo
        index = neighborCellInfo.index 
        neighborCell = cellGrid.cells[index[1], index[2], index[3]]
        append!(neighbors, neighborCell.atoms)
    end
    return neighbors
end




function GetTargetsFromNeighbor(atom::Atom, gridCell::GridCell, filterIndexes::Vector{Int64}, simulator::Simulator)
    cellGrid = simulator.cellGrid
    box = simulator.box
    targets = Vector{Atom}()
    infiniteFlag = true
    candidateTargets = Vector{Atom}()
    for (_, neighborCellInfo) in gridCell.neighborCellsInfo
        index = neighborCellInfo.index
        neighborCell = cellGrid.cells[index[1], index[2], index[3]]
        if neighborCell.isExplored
            continue
        end
        neighborCell.isExplored = true
        push!(simulator.exploredCells, neighborCell)
        infiniteFlag = false
        for neighborAtom in neighborCell.atoms
            if neighborAtom.index == atom.index || neighborAtom.index in filterIndexes    
                continue
            end
            Pertubation!(neighborAtom, simulator)
            if ComputeVDistance(atom, neighborAtom, neighborCellInfo.cross, box) > 0 
                ComputeP!(atom, neighborAtom, neighborCellInfo.cross, box)
                if neighborAtom.pValue[atom.index] >= simulator.parameters.pMax
                    DeleteP!(neighborAtom, atom.index)
                    if neighborAtom.latticePointIndex != -1
                        SetCoordinate!(neighborAtom, simulator.latticePoints[neighborAtom.latticePointIndex].coordinate)
                    end
                    continue
                end
                push!(candidateTargets, neighborAtom)
            end
        end
    end
    if infiniteFlag
        Log("Infinitely fly atom in the $(simulator.nCascade)th irradiation:\n$(atom)\n")
    end
    if isempty(candidateTargets)
        return (targets, infiniteFlag)
    end
    # Find target with minimum pL value using Julia's built-in findmin
    _, minIdx = findmin(t -> t.pL[atom.index], candidateTargets)
    nearestTarget = candidateTargets[minIdx]    
    push!(targets, nearestTarget)
    for candidateTarget in candidateTargets
        if candidateTarget.index == nearestTarget.index
            continue
        end
        matchFlag = true
        for target in targets            
            if !SimultaneousCriteria(atom, candidateTarget, target, simulator)
                matchFlag = false
                DeleteP!(candidateTarget, atom.index)
                if candidateTarget.latticePointIndex != -1
                    SetCoordinate!(candidateTarget, simulator.latticePoints[candidateTarget.latticePointIndex].coordinate)
                end
                break
            end
        end
        if matchFlag
            push!(targets, candidateTarget)
        end
    end    
    return (targets, infiniteFlag)
end


function DeleteP!(atom_t::Atom, atom_pIndex::Int64)
    delete!(atom_t.pValue, atom_pIndex)
    delete!(atom_t.pPoint, atom_pIndex)
    delete!(atom_t.pVector, atom_pIndex)
    delete!(atom_t.pL, atom_pIndex)
end

function EmptyP!(atom_t::Atom)
    empty!(atom_t.pValue)
    empty!(atom_t.pPoint)
    empty!(atom_t.pVector)
    empty!(atom_t.pL)
end


function SimultaneousCriteria(atom::Atom, neighborAtom::Atom, addedTarget::Atom, simulator::Simulator)
    flagP = abs(neighborAtom.pValue[atom.index] - addedTarget.pValue[atom.index]) <= simulator.parameters.pMax
    flagQ = abs(neighborAtom.pL[atom.index] - addedTarget.pL[atom.index]) <= simulator.constantsByType.qMax[[neighborAtom.type, addedTarget.type]]
    return flagP && flagQ
end


function ChangeCell!(atom::Atom, nextCellIndex::Vector{Int64}, simulator::Simulator)
    originalCell = simulator.cellGrid.cells[atom.cellIndex[1], atom.cellIndex[2], atom.cellIndex[3]]
    delete!(originalCell, atom, simulator)
    nextCell = simulator.cellGrid.cells[nextCellIndex[1], nextCellIndex[2], nextCellIndex[3]]
    push!(nextCell, atom, simulator)
end


function SetVelocityDirection!(atom::Atom, velocity::Vector{Float64})
    atom.velocityDirection = velocity / norm(velocity)
end

function SetEnergy!(atom::Atom, energy::Float64)
    if energy < 0.0
        atom.energy = 0.0
    else
        atom.energy = energy
    end
end

function GetNeighborVacancy(atom::Atom, simulator::Simulator)
    cells = simulator.cellGrid.cells        
    cell = cells[atom.cellIndex[1], atom.cellIndex[2], atom.cellIndex[3]]
    nearestVacancyDistance_squared = Inf
    nearestVacancyIndex = -1
    for (_, neighborCellInfo) in cell.neighborCellsInfo
        index = neighborCellInfo.index
        cross = neighborCellInfo.cross
        neighborCell = cells[index[1], index[2], index[3]]
        for latticePoint in neighborCell.latticePoints
            if latticePoint.atomIndex == -1
                dr2 = ComputeDistance_squared(atom.coordinate, latticePoint.coordinate, cross, simulator.box)
                if dr2 < simulator.parameters.vacancyRecoverDistance_squared && dr2 < nearestVacancyDistance_squared
                    nearestVacancyDistance_squared = dr2
                    nearestVacancyIndex = latticePoint.index
                end
            end 
        end
    end
    return nearestVacancyIndex
end


function Stop!(atom::Atom, simulator::Simulator)
    #SetEnergy!(atom, 0.0)
    Recover!(atom, simulator)
end


function Recover!(atom::Atom, simulator::Simulator)
    nearestVacancyIndex = GetNeighborVacancy(atom, simulator)
    if nearestVacancyIndex != -1
        SetOnLatticePoint!(atom, simulator.latticePoints[nearestVacancyIndex], simulator)
    end
end


function SetOnLatticePoint!(atom::Atom, latticePoint::LatticePoint, simulator::Simulator; isUpdateEnv::Bool = true)
    latticePoint.atomIndex = atom.index
    atom.latticePointIndex = latticePoint.index
    SetCoordinate!(atom, latticePoint.coordinate)
    if atom.isAlive && atom.cellIndex != latticePoint.cellIndex
        ChangeCell!(atom, latticePoint.cellIndex, simulator)
    elseif !atom.isAlive
        atom.isAlive = true
        nextCell = simulator.cellGrid.cells[latticePoint.cellIndex[1], latticePoint.cellIndex[2], latticePoint.cellIndex[3]]
        push!(nextCell, atom, simulator)
    end 
    if simulator.parameters.isKMC && isUpdateEnv
        latticePointIndexs = Set([latticePoint.environment;latticePoint.index])
        UpdateEvents!(latticePointIndexs, simulator)
    end
end


function AddToStore!(atom::Atom, simulator::Simulator)
    if simulator.isStore && atom.index <= simulator.numberOfAtomsWhenStored
        push!(simulator.displacedAtoms, atom.index)
    end 
end

function DeleteFromStore!(atom::Atom, simulator::Simulator)
    if simulator.isStore && atom.index <= simulator.numberOfAtomsWhenStored
        deleteat!(simulator.displacedAtoms, findfirst(==(atom.index), simulator.displacedAtoms))
    end
end 

function Restore!(simulator::Simulator)
    for atom in simulator.atoms[simulator.numberOfAtomsWhenStored+1:end]
        # Delete ions remained in the system from their cells.
        # Ions in simulator.atoms will be deleted latter by setting simulator.atoms = simulator.atoms[1:maxAtomID]. 
        if atom.isAlive
            delete!(simulator.cellGrid.cells[atom.cellIndex[1], atom.cellIndex[2], atom.cellIndex[3]], atom, simulator)
        end
    end
    for index in Set(simulator.displacedAtoms)
        atom = simulator.atoms[index]
        if atom.index == atom.latticePointIndex
            continue
        end
        latticePoint = simulator.latticePoints[atom.index]
        SetOnLatticePoint!(atom, latticePoint, simulator)
    end
    maxAtomID = simulator.numberOfAtomsWhenStored
    simulator.atoms = simulator.atoms[1:maxAtomID]
    simulator.maxAtomID = maxAtomID
    simulator.numberOfAtoms  = maxAtomID
    simulator.nCollisionEvent = 0
    empty!(simulator.displacedAtoms)
end


function Save!(simulator::Simulator)
    for atom in simulator.atoms
        if atom.latticePointIndex == -1
            error("Atom $(atom.index) is not on lattice when stored.")
        end
    end
    simulator.isStore = true
    simulator.numberOfAtomsWhenStored = simulator.numberOfAtoms
end 


function GetEnvironmentLatticePoints(latticePoint::LatticePoint, simulator::Simulator)
    cellIndex = latticePoint.cellIndex
    theCell = simulator.cellGrid.cells[cellIndex[1], cellIndex[2], cellIndex[3]]
    cells = simulator.cellGrid.cells
    cut_squared = simulator.environmentCut^2
    box = simulator.box
    environmentLatticePointsIndex = Vector{Int64}()
    dVectors = Vector{Vector{Float64}}()
    for (_, neighborCellInfo) in theCell.neighborCellsInfo
        index, cross = neighborCellInfo.index, neighborCellInfo.cross
        cell = cells[index[1], index[2], index[3]]
        latticePoints = cell.latticePoints
        for neighborLatticePoint in latticePoints
            neighborLatticePointIndex = neighborLatticePoint.index
            if ComputeDistance_squared(latticePoint.coordinate, neighborLatticePoint.coordinate, cross, box) <= cut_squared && neighborLatticePointIndex != latticePoint.index
                push!(dVectors, VectorDifference(latticePoint.coordinate, neighborLatticePoint.coordinate, cross, box))
                push!(environmentLatticePointsIndex, neighborLatticePointIndex)
            end
        end
    end
    # sort indexes by x then y then z of dVectors
    sorted_indices = sortperm(dVectors, by = v -> (v[1], v[2], v[3]))
    environmentLatticePointsIndex = environmentLatticePointsIndex[sorted_indices]
    
    return environmentLatticePointsIndex
end


function InitLatticePointEnvronment(simulator::Simulator)
    if simulator.parameters.DTEMode != 1 && simulator.parameters.DTEMode != 4
        println("Initializing lattice point environment...")
        for latticePoint in simulator.latticePoints
            latticePoint.environment = GetEnvironmentLatticePoints(latticePoint, simulator)
        end
    end
    # simulator.environmentLength should be get from the DTEDict.
end


function GetEnvironmentIndex(latticePoint::LatticePoint, simulator::Simulator)
    environment = latticePoint.environment
    latticePoints = simulator.latticePoints
    index = 0
    
    for i in 1:length(environment)
        if latticePoints[environment[i]].atomIndex != -1
            index += 2^(i-1)
        end
    end
    return index + 1
end 

function GaussianDeltaX(sigma::Float64)
    u1 = 0.0
    r2 = 0.0
    while true
        u1 = rand()*2-1
        u2 = rand()*2-1
        r2 = u1 * u1 + u2 * u2
        if r2 <= 1 && r2 > 0
            break
        end
    end
    return u1 * sqrt(-2 * log(r2) / r2) * sigma
end

function Pertubation!(atom::Atom, simulator::Simulator)
    if simulator.parameters.temperature > 0 
        for d in 1:3
            atom.coordinate[d] += GaussianDeltaX(simulator.constantsByType.sigma[atom.type])
        end
    end
end

function SetCoordinate!(atom::Atom, coordinate::Vector{Float64})
    atom.coordinate .= coordinate
end 



function TemperatureToSigma(T::Float64, ΘD::Float64, M_u::Float64)
    ħ  = 1.054_571_817e-34   # J·s
    kB = 1.380_649e-23       # J/K
    Å  = 1e-10               # meter
    M = M_u * 1.660_539_066e-27             
    a  = ħ^2 / (2*M*kB*ΘD)
    y  = ΘD / T
    integ = quadgk(x -> x/(exp(x)-1), 0, y)[1]
    u2 = a * (0.25 + (T/ΘD)^2 * integ)       # m²
    return sqrt(u2) / Å                     # Å
end




