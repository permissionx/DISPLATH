using StaticArrays

function Box(Vectors::Matrix{Float64})
    log_info("Box created: $(round(Vectors[1,1]; digits=2)) × $(round(Vectors[2,2]; digits=2)) × $(round(Vectors[3,3]; digits=2)) Å")
    return Box(Vectors, inv(Vectors'), true)
end 

function CreateBoxByPrimaryVectors(primaryVectors::Matrix{Float64}, sizes::Vector{Int64})
    vectors = primaryVectors .* sizes
    return Box(vectors)
end 


function Atom(type::Int64, coordinate::Vector{Float64}, parameters::Parameters)
    index = 0
    isAlive = true
    cellIndex = (0,0,0)
    velocityDirection = SVector{3,Float64}(0.0, 0.0, 0.0)  
    energy = 0.0
    radius, mass, Z, dte, bde, _, _ = TypeToProperties(type, parameters.typeDict)
    #numberOfEmptyCells = 0
    emptyPath = 0.0
    pValue = 0.0
    pVector = SVector{3,Float64}(0.0, 0.0, 0.0)  
    pPoint = SVector{3,Float64}(0.0, 0.0, 0.0)   
    lastTargets = Vector{Int64}()
    pL = 0.0
    pAtomIndex = -1 # temperory 
    pDirection = Float64[0.0,0.0,0.0] # temperory 
    latticePointIndex = -1
    frequency = 0.0
    frequencies = Vector{Float64}()
    finalLatticePointEnvIndexs = Vector{Int64}()
    eventIndex = -1
    isNewlyLoaded = false
    lattcieCoordinate = SVector{3,Float64}(coordinate[1], coordinate[2], coordinate[3])  
    indexInCell = 0
    return Atom(index, isAlive, type, coordinate[:], cellIndex, 
                radius, mass, velocityDirection, energy, Z, 
                dte, bde, emptyPath, #numberOfEmptyCells,
                pValue, pVector, pPoint, pL, pAtomIndex, pDirection, lastTargets, # temperory 
                latticePointIndex,
                frequency, frequencies, finalLatticePointEnvIndexs, eventIndex, 
                isNewlyLoaded, lattcieCoordinate, indexInCell)
end


function TypeToProperties(type::Int64, typeDict::Dict{Int64, Element})
    if haskey(typeDict, type)
        element = typeDict[type]
        return element.radius, element.mass, element.Z, element.dte, element.bde, element.alpha, element.beta 
    else
        error("Unknown atom type: $type")
    end 
end 


function SetCellNeighborInfo!(cell::Cell, grid::Grid)
    # Direct triple loop implementation - much faster than recursion
    for delta_x in [-1, 0, 1]
        for delta_y in [-1, 0, 1]
            for delta_z in [-1, 0, 1]
                neighborKeys = (Int8(delta_x), Int8(delta_y), Int8(delta_z))  
                neighborIndex = [0, 0, 0]  
                neighborCross = [Int8(0), Int8(0), Int8(0)]  
                # Calculate neighbor cell index and cross flags for each dimension
                for d in 1:3
                    delta = neighborKeys[d]
                    index = cell.index[d] + delta
                    cross = Int8(0)
                    if index < 1
                        index += grid.sizes[d]
                        cross = Int8(-1)
                    elseif index > grid.sizes[d]
                        index -= grid.sizes[d]
                        cross = Int8(1)
                    end
                    neighborIndex[d] = index
                    neighborCross[d] = cross
                end
                neighborIndex_tuple = (neighborIndex[1], neighborIndex[2], neighborIndex[3])
                neighborCross_tuple = (neighborCross[1], neighborCross[2], neighborCross[3])
                neighborCellInfo = NeighborCellInfo(neighborIndex_tuple, neighborCross_tuple)
                idx = (delta_x + 2, delta_y + 2, delta_z + 2)
                cell.neighborCellsInfo[idx...] = neighborCellInfo
            end
        end
    end
end


function CreateGrid(box::Box, inputVectors::Matrix{Float64})
    if !box.isOrthogonal
        error("The box is not orthogonal, please use the orthogonal box.")
    end
    sizes = Vector{Int64}(undef, 3)
    vectors = Matrix{Float64}(undef, 3, 3)
    for d in 1:3
        sizes[d] = Int64(floor(box.vectors[d,d] / inputVectors[d,d]))
        vectors[d,d] = box.vectors[d,d] / sizes[d]
    end
    log_info("Cell grid: $(sizes[1]) × $(sizes[2]) × $(sizes[3]) = $(sizes[1]*sizes[2]*sizes[3]) cells")
    log_info("Cell size: $(round(vectors[1,1]; digits=2)) × $(round(vectors[2,2]; digits=2)) × $(round(vectors[3,3]; digits=2)) Å")
    if ! IS_DYNAMIC_LOAD
        cells = Array{Cell, 3}(undef, sizes[1], sizes[2], sizes[3])
        @showprogress desc="Creating cells: " for x in 1:sizes[1]
            for y in 1:sizes[2]
                for z in 1:sizes[3]
                    cells[x, y, z] = CreateCell((x, y, z), vectors)
                end
            end    
        end
        cellVolume = vectors[1,1] * vectors[2,2] * vectors[3,3]
        grid = Grid(cells, vectors, sizes, cellVolume) 
        @showprogress desc="Pushing cell neighbors: " for cell in grid.cells
            SetCellNeighborInfo!(cell, grid)
        end
    else
        cells = Dict{Tuple{Int64, Int64, Int64}, Cell}()    
        cellVolume = vectors[1,1] * vectors[2,2] * vectors[3,3]
        grid = Grid(cells, vectors, sizes, cellVolume) 
    end
    log_success("Cell grid created")
    log_separator()
    return grid
end


function _GetCellDense(grid::Grid, cellIndex::Tuple{Int64, Int64, Int64})
    return grid.cells[cellIndex...]
end


function CreateCell(cellIndex::Tuple{Int64, Int64, Int64}, vectors::Matrix{Float64})
    x, y, z = cellIndex
    ranges = Matrix{Float64}(undef, 3, 2)
    ranges[1,1] = (x-1) * vectors[1,1]
    ranges[1,2] = x * vectors[1,1]
    ranges[2,1] = (y-1) * vectors[2,2]
    ranges[2,2] = y * vectors[2,2]
    ranges[3,1] = (z-1) * vectors[3,3]
    ranges[3,2] = z * vectors[3,3]  
    cell = Cell(cellIndex, Vector{Atom}(), Vector{LatticePoint}(), 
                            ranges, 
                            Array{NeighborCellInfo, 3}(undef, 3, 3, 3), false, 0.0)
    return cell
end


function _GetCellDict!(grid::Grid, cellIndex::Tuple{Int64, Int64, Int64})
    return get!(grid.cells, cellIndex) do
        CreateCell(cellIndex, grid.vectors)
    end
end 


function GetCell(grid::Grid, cellIndex::Tuple{Int64, Int64, Int64})
    if !IS_DYNAMIC_LOAD
        return _GetCellDense(grid, cellIndex)
    else
        return _GetCellDict!(grid, cellIndex)
    end
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
    log_info("")
    log_info("Vibration σ for each type:")
    for p in types
        radius_p, mass_p, Z_p, _, _, α_p, β_p = TypeToProperties(p, typeDict)
        for t in types
            radius_t, _, Z_t, _, _, _, _ = TypeToProperties(t, typeDict)
            V_upterm[[p,t]] = BCA.ConstantFunctions.V_upterm(Z_p, Z_t)
            a_U[[p,t]] = BCA.ConstantFunctions.a_U(Z_p, Z_t)
            S_e_upTerm[[p,t]] = BCA.ConstantFunctions.S_e_upTerm(p, Z_p, Z_t, mass_p, α_p)
            x_nl[[p,t]] = BCA.ConstantFunctions.x_nl(p, Z_p, Z_t, β_p)
            a[[p,t]] = BCA.ConstantFunctions.a(Z_p, Z_t)
            Q_nl[[p,t]] = BCA.ConstantFunctions.Q_nl(Z_p, Z_t, parameters.pMax)
            Q_loc[[p,t]] = BCA.ConstantFunctions.Q_loc(Z_p, Z_t)
            qMax[[p,t]] = radius_p + radius_t
        end
        E_m[p] = BCA.ConstantFunctions.E_m(Z_p, mass_p)
        sigma[p] = TemperatureToSigma(parameters.temperature, parameters.DebyeTemperature, mass_p)
        log_info("  Type $(p): σ = $(round(sigma[p]; digits=3)) Å")
    end
    return ConstantsByType(V_upterm, a_U, E_m, S_e_upTerm, S_e_downTerm, x_nl, a, Q_nl, Q_loc, qMax, sigma)
end


function InitθτFunctions(parameters::Parameters, constantsByType::ConstantsByType)
    typeDict = parameters.typeDict
    θFunctions = Dict{Vector{Int64}, Function}()
    τFunctions = Dict{Vector{Int64}, Function}()
    log_separator()
    log_info("Loading θ and τ functions...")
    for type_p in keys(typeDict)
        for type_t in keys(typeDict)
            mass_p = typeDict[type_p].mass
            mass_t = typeDict[type_t].mass
            θInterpolation, τInterpolation = θτFunctions(mass_p, mass_t, type_p, type_t, constantsByType, parameters)
            θFunctions[[type_p, type_t]] = (E_p, p) -> θInterpolation(E_p, p)
            τFunctions[[type_p, type_t]] = (E_p, p) -> τInterpolation(E_p, p)
            log_debug("  $(parameters.typeDict[type_p].name) → $(parameters.typeDict[type_t].name) loaded")
        end
    end
    log_success("All θ and τ functions initialized")
    log_separator()
    return θFunctions, τFunctions
end


function θτFunctions(mass_p::Float64, mass_t::Float64, type_p::Int64, type_t::Int64, constantsByType::ConstantsByType, parameters::Parameters)
    E_p_axis = Float64[]
    p_axis = Float64[]
    θMatrix = Matrix{Float64}(undef, 0, 0)
    τMatrix = Matrix{Float64}(undef, 0, 0)
    try
        E_p_axis, p_axis, θMatrix, τMatrix = LoadθτData(type_p, type_t, parameters)
    catch
        EPowerRange = parameters.EPowerRange    
        pPowerRange = parameters.pPowerRange
        nE = length(EPowerRange)
        np = length(pPowerRange)
        θMatrix = Array{Float64, 2}(undef, nE, np)
        τMatrix = Array{Float64, 2}(undef, nE, np)
        N = length(EPowerRange)
        @showprogress @threads for i in 1:N
            E_p_power = EPowerRange[i]
            E_p = 10.0^E_p_power
            for (j, p_power) in enumerate(pPowerRange)
                p = 10.0^p_power
                θ, τ = BCA.θτ(E_p, mass_p, mass_t, type_p, type_t, p, constantsByType)
                θMatrix[i, j] = θ
                τMatrix[i, j] = τ
            end
        end
        E_p_axis = collect(EPowerRange)
        p_axis = collect(pPowerRange)    
        SaveθτData(type_p, type_t, θMatrix, τMatrix, E_p_axis, p_axis, parameters)
    end
    # interpolate
    θFunction = interpolate((E_p_axis, p_axis), θMatrix, Gridded(Linear()))
    τFunction = interpolate((E_p_axis, p_axis), τMatrix, Gridded(Linear()))
    return θFunction, τFunction
end


function Simulator(box::Box, atoms::Vector{Atom}, inputGridVectors::Matrix{Float64}, parameters::Parameters)
    log_section("Initializing Simulator.")
    simulator = Simulator(box, inputGridVectors, parameters)
    if !IS_DYNAMIC_LOAD
        LoadAtoms!(simulator, atoms)
    end
    log_success("Simulator initialized.")
    return simulator 
end

function LoadAtoms!(simulator::Simulator, atoms::Vector{Atom})
    for atom in atoms
        push!(simulator, atom) 
        latticePoint = LatticePoint(atom)
        push!(simulator, latticePoint)
    end
    for cell in simulator.grid.cells
        cell.atomicDensity = length(cell.atoms) / simulator.grid.cellVolume
    end 
    for atom in simulator.atoms
        Pertubation!(atom, simulator)
    end
    InitLatticePointEnvronment(simulator)
    log_success("$(simulator.numberOfAtoms) atoms loaded.")
end

function CreateAtomsByPrimaryVectors(parameters::Parameters)
    primaryVectors = parameters.primaryVectors
    latticeRanges = parameters.latticeRanges
    basis = parameters.basis
    basisTypes = parameters.basisTypes
    atomNumber = (latticeRanges[1,2] - latticeRanges[1,1]) * (latticeRanges[2,2] - latticeRanges[2,1]) * (latticeRanges[3,2] - latticeRanges[3,1]) * length(basisTypes)
    atoms = Vector{Atom}(undef, atomNumber)
    n = 1
    @showprogress desc="Creating atoms ($(atomNumber)): " for x in latticeRanges[1,1]:latticeRanges[1,2]-1
        for y in latticeRanges[2,1]:latticeRanges[2,2]-1    
            for z in latticeRanges[3,1]:latticeRanges[3,2]-1
                for i in 1:length(basisTypes)
                    reducedCoordinate = Float64[x,y,z] + basis[i, :]
                    coordinate = primaryVectors' * reducedCoordinate
                    atoms[n] = Atom(basisTypes[i], coordinate, parameters)
                    n += 1
                end
            end
        end
    end
    return atoms
end


function Simulator(boxVectors::Matrix{Float64}, inputGridVectors::Matrix{Float64}, parameters::Parameters)
    box = Box(boxVectors)
    if !IS_DYNAMIC_LOAD
        atoms = CreateAtomsByPrimaryVectors(parameters)
    end
    simulator = Simulator(box, atoms, inputGridVectors, parameters)
    return simulator    
end 


function Simulator(boxSizes::Vector{Int64}, inputGridVectors::Matrix{Float64}, parameters::Parameters)
    box = CreateBoxByPrimaryVectors(parameters.primaryVectors, boxSizes)
    if !IS_DYNAMIC_LOAD
        atoms = CreateAtomsByPrimaryVectors(parameters)
    end
    simulator = Simulator(box, atoms, inputGridVectors, parameters)
    return simulator    
end


function Parameters(pMax::Float64, vacancyRecoverDistance::Float64, typeDict::Dict{Int64, Element}; kwargs...)
    primaryVectors = [1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0]
    latticeRanges = [0 1; 0 1; 0 1]
    basis = [0.0 0.0 0.0]
    basisTypes = [1]  
    parameters = Parameters(primaryVectors, latticeRanges, basisTypes, basis, pMax, vacancyRecoverDistance, typeDict; kwargs...)
    return parameters
end

function LoadAtomsAndBoxFromDataFile(fileName::String; replicate::Vector{Int64} = [1,1,1])
    xlo, xhi, ylo, yhi, zlo, zhi, types, xs, ys, zs = ReadDate(fileName, replicate)
    box = Box([xhi-xlo 0.0 0.0; 0.0 yhi-ylo 0.0; 0.0 0.0 zhi-zlo])
    atoms = Vector{Atom}(undef, length(types))
    for (n,(type, x, y, z)) in enumerate(zip(types, xs, ys, zs))
        atoms[n] = Atom(type, [x, y, z], parameters)
    end
    return box, atoms
end

function Simulator(fileName::String, inputGridVectors::Matrix{Float64}, parameters::Parameters; replicate::Vector{Int64} = [1,1,1])
    if IS_DYNAMIC_LOAD
        error("Simulator from date file is not supported in dynamic load mode.")
    end 
    box, atoms = LoadAtomsAndBoxFromDataFile(fileName; replicate=replicate)
    simulator = Simulator(box, atoms, inputGridVectors, parameters)
    return simulator
end


function LatticePoint(atom::Atom)
    environment = Vector{Int64}() 
    return LatticePoint(copy(atom.index), copy(atom.type), 
                        copy(atom.coordinate), atom.cellIndex, environment,
                        atom.index)
end


function WhichCell(coordinate::Vector{Float64}, grid::Grid)
    cellIndex = Vector{Int64}(undef, 3)
    for d in 1:3
        cellIndex[d] = Int64(floor(coordinate[d] / grid.vectors[d,d])) + 1
        if cellIndex[d] < 1 
            cellIndex[d] = 1
        elseif cellIndex[d] > grid.sizes[d]
            cellIndex[d] = grid.sizes[d]
        end
    end
    return (cellIndex[1], cellIndex[2], cellIndex[3])
end


function push!(simulator::Simulator, atom::Atom)
    atom.index = simulator.maxAtomID + 1
    simulator.maxAtomID += 1
    push!(simulator.atoms, atom)
    simulator.numberOfAtoms += 1
    cellIndex = WhichCell(atom.coordinate, simulator.grid)
    atom.cellIndex = cellIndex
    push!(GetCell(simulator.grid, cellIndex).atoms, atom)
end 


function push!(simulator::Simulator, latticePoint::LatticePoint)
    push!(simulator.latticePoints, latticePoint)
    push!(GetCell(simulator.grid, latticePoint.cellIndex).latticePoints, latticePoint)
    simulator.atoms[latticePoint.atomIndex].latticePointIndex = latticePoint.index
end 

function push!(cell::Cell, atom::Atom, simulator::Simulator)
    atom.cellIndex = cell.index
    push!(cell.atoms, atom)
    if !IS_DYNAMIC_LOAD
        cell.atomicDensity = length(cell.atoms) / simulator.grid.cellVolume
    end
end 


function delete!(simulator::Simulator, atom::Atom)
    originalCell = GetCell(simulator.grid, atom.cellIndex)
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
    vacancy = Atom(latticePoint.type, latticePoint.coordinate, simulator.parameters)
    vacancy.index = latticePoint.index
    push!(simulator.vacancies, vacancy)
    push!(GetCell(simulator.grid, latticePoint.cellIndex).vacancies, vacancy)
    if isUpdateEnv && simulator.parameters.isKMC
        DeleteAtomEvents!(simulator, atom)
        UpdateEvents!(Set(latticePoint.environment), simulator)
    end
end


function delete!(cell::Cell, atom::Atom, simulator::Simulator)
    if !atom.isAlive
        error("Atom $(atom.index) is not alive when deleting")
    end
    #@show atom.index, atom.cellIndex, simulator.nCascade, simulator.nCollisionEvent, atom.coordinate
    deleteat!(cell.atoms, findfirst(a -> a.index == atom.index, cell.atoms))
    atom.cellIndex = (-1, -1, -1)
    if !IS_DYNAMIC_LOAD
        cell.atomicDensity = length(cell.atoms) / simulator.grid.cellVolume  
    end
end


function DisplaceAtom!(atom::Atom, newPosition::SVector{3,Float64}, simulator::Simulator)
    pos = if newPosition isa SVector
        [newPosition[1], newPosition[2], newPosition[3]]
    else
        copy(newPosition)
    end
    
    for d in 1:3
        # need to adappt non-periodic condition
        if pos[d] < 0
            if simulator.parameters.periodic[d] == false
                pos[d] = 0.01
            else
                pos[d] += simulator.box.vectors[d,d]
            end
        elseif pos[d] >= simulator.box.vectors[d,d]
            if simulator.parameters.periodic[d] == false
                pos[d] = simulator.box.vectors[d,d] - 0.01
            else
                pos[d] -= simulator.box.vectors[d,d]
            end
        end
    end
    
    SetCoordinate!(atom, pos)
    cellIndex = WhichCell(atom.coordinate, simulator.grid)

    if cellIndex != atom.cellIndex
        ChangeCell!(atom, cellIndex, simulator)
    end
end

function DisplaceAtom!(atom::Atom, newPosition::Vector{Float64}, simulator::Simulator)
    DisplaceAtom!(atom, SVector{3,Float64}(newPosition[1], newPosition[2], newPosition[3]), simulator)
end


function ComputeDistance_squared(coordinate1::Vector{Float64}, coordinate2::Vector{Float64}, crossFlag::NTuple{3, Int8}, box::Box)
    dv = VectorDifference(coordinate1, coordinate2, crossFlag, box)
    distance_squared = dv[1]* dv[1] + dv[2]*dv[2] + dv[3]  * dv[3]
    return distance_squared
end


function ComputeDistance(coordinate1::Vector{Float64}, coordinate2::Vector{Float64}, crossFlag::NTuple{3, Int8}, box::Box)
    return sqrt(ComputeDistance_squared(coordinate1, coordinate2, crossFlag, box))
end


function ComputeVDistance(atom_p::Atom, atom_t::Atom, crossFlag::NTuple{3, Int8}, box::Box)
    # v for atom_p
    dv = VectorDifference(atom_p.coordinate, atom_t.coordinate, crossFlag, box)
    return dv' * atom_p.velocityDirection
end


function VectorDifference(v1::Vector{Float64}, v2::Vector{Float64}, crossFlag::NTuple{3, Int8}, box::Box)
    if crossFlag == (Int8(0), Int8(0), Int8(0))
        return v2 - v1
    end 
    return SVector{3,Float64}(
        v2[1] - v1[1] + crossFlag[1] * box.vectors[1,1],
        v2[2] - v1[2] + crossFlag[2] * box.vectors[2,2],
        v2[3] - v1[3] + crossFlag[3] * box.vectors[3,3]
    )
end


function ComputeP!(atom_p::Atom, atom_t::Atom, crossFlag::NTuple{3, Int8}, box::Box)
    dv = VectorDifference(atom_p.coordinate, atom_t.coordinate, crossFlag, box)
    t = dot(dv, atom_p.velocityDirection)
    atom_t.pL = t
    if 1 in crossFlag || -1 in crossFlag
        pPoint_calc = Vector{Float64}(atom_p.coordinate + t * atom_p.velocityDirection)
        for d in 1:3
            if crossFlag[d] != 0
                pPoint_calc[d] -= crossFlag[d] * box.vectors[d,d]
            end
        end
    else
        pPoint_calc = atom_p.coordinate + t * atom_p.velocityDirection
    end
    atom_t.pPoint = SVector{3,Float64}(pPoint_calc[1], pPoint_calc[2], pPoint_calc[3])
    pVector_calc = atom_t.pPoint - atom_t.coordinate
    atom_t.pVector = SVector{3,Float64}(pVector_calc[1], pVector_calc[2], pVector_calc[3])
    p = norm(atom_t.pVector)
    atom_t.pValue = p
    # need to check periodic condition
    return p
end




function ChangeCell!(atom::Atom, nextCellIndex::Tuple{Int64, Int64, Int64}, simulator::Simulator)
    originalCell = GetCell(simulator.grid, atom.cellIndex)
    delete!(originalCell, atom, simulator)
    nextCell = GetCell(simulator.grid, nextCellIndex)
    push!(nextCell, atom, simulator)
end


function SetVelocityDirection!(atom::Atom, velocity::SVector{3,Float64})
    n = norm(velocity)
    if isnan(n) || n == Inf || n == 0.0
        atom.velocityDirection = SVector{3,Float64}(0.0, 0.0, 0.0)
    else
        normalized_velocity = velocity / n
        atom.velocityDirection = SVector{3,Float64}(normalized_velocity[1], normalized_velocity[2], normalized_velocity[3])
    end
end

function SetVelocityDirection!(atom::Atom, velocity::Vector{Float64})
    SetVelocityDirection!(atom, SVector{3,Float64}(velocity[1], velocity[2], velocity[3]))
end


function SetEnergy!(atom::Atom, energy::Float64)
    if energy < 0.0
        atom.energy = 0.0
    else
        atom.energy = energy
    end
end

function GetNeighborVacancy(atom::Atom, simulator::Simulator)
    grid = simulator.grid    
    cell = GetCell(grid, atom.cellIndex)
    nearestVacancyDistance_squared = Inf
    nearestVacancyIndex = -1
    for neighborCellInfo in cell.neighborCellsInfo
        index, cross = neighborCellInfo.index, neighborCellInfo.cross
        neighborCell = GetCell(grid, index)
        for vacancy in neighborCell.vacancies
            dr2 = ComputeDistance_squared(atom.coordinate, vacancy.coordinate, cross, simulator.box)
            if dr2 < simulator.parameters.vacancyRecoverDistance_squared && dr2 < nearestVacancyDistance_squared
                nearestVacancyDistance_squared = dr2
                nearestVacancyIndex = vacancy.index
            end
        end
    end
    return nearestVacancyIndex
end


function Stop!(atom::Atom, simulator::Simulator)
    SetEnergy!(atom, 0.0)
    SetVelocityDirection!(atom, SVector{3,Float64}([0.0, 0.0, 0.0]))
    Recover!(atom, simulator)
end


function Recover!(atom::Atom, simulator::Simulator)
    nearestVacancyIndex = GetNeighborVacancy(atom, simulator)
    if nearestVacancyIndex != -1
        SetOnLatticePoint!(atom, simulator.latticePoints[nearestVacancyIndex], simulator)
        deleteat!(simulator.vacancies, findfirst(v -> v.index == nearestVacancyIndex, simulator.vacancies))
        cell = GetCell(simulator.grid, atom.cellIndex)
        deleteat!(cell.vacancies, findfirst(v -> v.index == nearestVacancyIndex, cell.vacancies))
    end
end


function SetOnLatticePoint!(atom::Atom, latticePoint::LatticePoint, simulator::Simulator; isUpdateEnv::Bool = true)
    SetEnergy!(atom, 0.0)
    SetVelocityDirection!(atom, SVector{3,Float64}([0.0, 0.0, 0.0]))
    latticePoint.atomIndex = atom.index
    atom.latticePointIndex = latticePoint.index
    SetCoordinate!(atom, latticePoint.coordinate)
    if atom.isAlive && atom.cellIndex != latticePoint.cellIndex
        ChangeCell!(atom, latticePoint.cellIndex, simulator)
    elseif !atom.isAlive
        atom.isAlive = true
        nextCell = GetCell(simulator.grid, latticePoint.cellIndex)
        push!(nextCell, atom, simulator)
    end 
    if simulator.parameters.isKMC && isUpdateEnv
        latticePointIndexs = Set([latticePoint.environment;latticePoint.index])
        UpdateEvents!(latticePointIndexs, simulator)
    end
    Pertubation!(atom, simulator)
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
    if ! IS_DYNAMIC_LOAD
        Restore_staticLoad!(simulator)
    else
        Restore_dynamicLoad!(simulator)
    end
end

function Restore_staticLoad!(simulator::Simulator)
    for atom in simulator.atoms[simulator.numberOfAtomsWhenStored+1:end]
        # Delete ions remained in the system from their cells.
        # Ions in simulator.atoms will be deleted latter by setting simulator.atoms = simulator.atoms[1:maxAtomID]. 
        if atom.isAlive
            delete!(GetCell(simulator.grid, atom.cellIndex), atom, simulator)
        end
    end
    latticePoints = simulator.latticePoints
    for vacancy in simulator.vacancies
        latticePoint = latticePoints[vacancy.index]
        cell = GetCell(simulator.grid, latticePoint.cellIndex)
        empty!(cell.vacancies)
    end
    empty!(simulator.vacancies)

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
    theCell = GetCell(simulator.grid, cellIndex)
    grid = simulator.grid
    cut_squared = simulator.environmentCut^2
    box = simulator.box
    environmentLatticePointsIndex = Vector{Int64}()
    dVectors = Vector{Vector{Float64}}()
    for neighborCellInfo in theCell.neighborCellsInfo
        index, cross = neighborCellInfo.index, neighborCellInfo.cross
        cell = GetCell(grid, index)
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
        log_info("🌐 Initializing lattice point environment...\n")
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
    return randn(THREAD_RNG[Threads.threadid()]) * sigma 
end


function Pertubation!(atom::Atom, simulator::Simulator)
    if simulator.parameters.isAmorphous 
        rng = THREAD_RNG[Threads.threadid()]
        atom.coordinate .= GetCell(simulator.grid, atom.cellIndex).ranges[:,1] .+ [rand(rng) * simulator.grid.vectors[d, d] for d in 1:3]
    else
        ah = simulator.parameters.amorphousHeight
        if atom.coordinate[3] > ah
            rng = THREAD_RNG[Threads.threadid()]
            cell = GetCell(simulator.grid, atom.cellIndex)
            atom.coordinate[1] = cell.ranges[1,1] + rand(rng) * simulator.grid.vectors[1, 1]
            atom.coordinate[2] = cell.ranges[2,1] + rand(rng) * simulator.grid.vectors[2, 2]
            base = cell.ranges[3,1] > ah ? cell.ranges[3,1] : ah
            latticeTop = simulator.parameters.primaryVectors[3,3] * simulator.parameters.latticeRanges[3,2]     
            top = cell.ranges[3,2] < latticeTop ? cell.ranges[3,2] : latticeTop
            atom.coordinate[3] = base + rand(rng) * (top - base)
        else
            if simulator.parameters.temperature > 0.0
                for d in 1:3
                    atom.coordinate[d] += GaussianDeltaX(simulator.constantsByType.sigma[atom.type])
                end
            end
        end
    end
end


function SetCoordinate!(atom::Atom, coordinate::Vector{Float64})
    atom.coordinate .= coordinate
end 


function TemperatureToSigma(T::Float64, θ_D::Float64, m_rel::Float64; atol=1e-10, rtol=1e-8)
    if T == 0.0
        log_debug("Temperature is 0 K")
        return 0
    end
    ħ   = 1.054_571_817e-34      # J·s
    kB  = 1.380_649_000e-23      # J/K
    amu = 1.660_539_066_60e-27   # kg

    M = m_rel * amu
    y_max = θ_D / T      

    # Integration of x/(e^x-1) dx
    integrand(x) = x / (exp(x) - 1)
    I, _ = quadgk(integrand, 0.0, y_max; atol, rtol)

    σ2 = 3 * ħ^2 / (M * kB * θ_D) * (0.25 + (T/θ_D)^2 * I)
    σ  = sqrt(σ2) * 1e10         # m → Å

    return σ
end
