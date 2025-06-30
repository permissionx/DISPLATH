function Box(Vectors::Matrix{Float64})
    # println("\n🚀 Initializing the simulator...\n")  # 删除多余的初始化输出
log_success("📦 Box created! Size: $(round(Vectors[1,1]; digits=3)) × $(round(Vectors[2,2]; digits=3)) × $(round(Vectors[3,3]; digits=3)) Å\n")
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
    velocityDirection = Float64[0.0,0.0,0.0]  # length： 0 or one
    energy = 0.0
    radius, mass, Z, dte, bde, _, _ = TypeToProperties(type, parameters.typeDict)
    numberOfEmptyCells = 0
    pValue = 0.0
    pVector = Vector{Float64}(undef, 3)
    pPoint = Vector{Float64}(undef, 3)  
    lastTargets = Vector{Int64}()
    pL = 0.0
    latticePointIndex = -1
    frequency = 0.0
    frequencies = Vector{Float64}()
    finalLatticePointEnvIndexs = Vector{Int64}()
    eventIndex = -1
    isNewlyLoaded = false
    lattcieCoordinate = Vector{Float64}(undef, 3)
    return Atom(index, isAlive, type, coordinate, cellIndex, 
                radius, mass, velocityDirection, energy, Z, 
                dte, bde, numberOfEmptyCells,
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


function SetCellNeighborInfo!(cell::Cell, grid::Grid)
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
                    index = cell.index[d] + delta
                    cross = 0
                    if index < 1
                        index += grid.sizes[d]
                        cross = -1
                    elseif index > grid.sizes[d]
                        index -= grid.sizes[d]
                        cross = 1
                    end
                    neighborIndex[d] = index
                    neighborCross[d] = cross
                end
                neighborIndex = Tuple(neighborIndex)
                neighborCross = Tuple(neighborCross)
                neighborCellInfo = NeighborCellInfo(neighborIndex, neighborCross)
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
    log_info("🧩 Cell grid: $(sizes[1]) × $(sizes[2]) × $(sizes[3]) = $(sizes[1]*sizes[2]*sizes[3]) cells, each size: $(round(vectors[1,1]; digits=3)) × $(round(vectors[2,2]; digits=3)) × $(round(vectors[3,3]; digits=3)) Å\n")
    log_info("----------------------------------------")
    if ! IS_DYNAMIC_LOAD
        cells = Array{Cell, 3}(undef, sizes[1], sizes[2], sizes[3])
        @showprogress desc="Creating cells: " for x in 1:sizes[1]
            for y in 1:sizes[2]
                for z in 1:sizes[3]
                    ranges = Matrix{Float64}(undef, 3, 2)
                    ranges[1,1] = (x-1) * vectors[1,1]
                    ranges[1,2] = x * vectors[1,1]
                    ranges[2,1] = (y-1) * vectors[2,2]
                    ranges[2,2] = y * vectors[2,2]
                    ranges[3,1] = (z-1) * vectors[3,3]
                    ranges[3,2] = z * vectors[3,3]  
                    cells[x, y, z] = Cell((x, y, z), Vector{Atom}(), Vector{LatticePoint}(), 
                                             ranges, Array{NeighborCellInfo, 3}(undef, 3, 3, 3), 
                                             false, 0.0)
                end
            end    
        end
        cellVolume = vectors[1,1] * vectors[2,2] * vectors[3,3]
        grid = Grid(cells, vectors, sizes, cellVolume) 
        @showprogress desc="Pushing cell neighbors: " for cell in grid.cells
            SetCellNeighborInfo!(cell, grid)
        end
    else
        cells = spzeros(Cell, sizes[1] * sizes[2] * sizes[3])
        cellVolume = vectors[1,1] * vectors[2,2] * vectors[3,3]
        grid = Grid(cells, vectors, sizes, cellVolume) 
    end
    log_success("✅ Cell grid created!\n----------------------------------------\n")
    return grid
end


function _GetCellDense(grid::Grid, cellIndex::Tuple{Int64, Int64, Int64})
    return grid.cells[cellIndex...]
end

function LinearIndex(cellIndex::Tuple{Int64, Int64, Int64}, sizes::Vector{Int64})
    x, y, z = cellIndex
    NX, NY = sizes[1], sizes[2]
    return (x-1) + (y-1)*NX + (z-1)*NX*NY + 1   # 1-based
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

function _GetCellSparse!(grid::Grid, cellIndex::Tuple{Int64, Int64, Int64})
    idx = LinearIndex(cellIndex, grid.sizes)
    if grid.cells[idx] === EMPTY_GRIDCELL
        cell = CreateCell(cellIndex, grid.vectors)
        grid.cells[idx] = cell
    end
    return grid.cells[idx]
end

function GetCell(grid::Grid, cellIndex::Tuple{Int64, Int64, Int64})
    if !IS_DYNAMIC_LOAD
        return _GetCellDense(grid, cellIndex)
    else
        return _GetCellSparse!(grid, cellIndex)
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
        log_info("  • Type $(p): σ = $(round(sigma[p]; digits=4)) Å")
    end
    log_info("\n----------------------------------------")
    return ConstantsByType(V_upterm, a_U, E_m, S_e_upTerm, S_e_downTerm, x_nl, a, Q_nl, Q_loc, qMax, sigma)
end


function InitθτFunctions(parameters::Parameters, constantsByType::ConstantsByType)
    typeDict = parameters.typeDict
    θFunctions = Dict{Vector{Int64}, Function}()
    τFunctions = Dict{Vector{Int64}, Function}()
    log_info("📥 θ and τ functions loading:")
    for type_p in keys(typeDict)
        for type_t in keys(typeDict)
            mass_p = typeDict[type_p].mass
            mass_t = typeDict[type_t].mass
            θInterpolation, τInterpolation = θτFunctions(mass_p, mass_t, type_p, type_t, constantsByType, parameters)
            θFunctions[[type_p, type_t]] = (E_p, p) -> θInterpolation(E_p, p)
            τFunctions[[type_p, type_t]] = (E_p, p) -> τInterpolation(E_p, p)
            log_success("  📄 $(parameters.typeDict[type_p].name) → $(parameters.typeDict[type_t].name) loaded. ✅")
        end
    end
    log_success("\n🧩 All θ and τ functions initialized! 🎉\n----------------------------------------\n")
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
        pRange = parameters.pRange
        nE = length(EPowerRange)
        np = length(pRange)
        θMatrix = Array{Float64, 2}(undef, nE, np)
        τMatrix = Array{Float64, 2}(undef, nE, np)
        N = length(EPowerRange)
        @showprogress @threads for i in 1:N
            E_p_power = EPowerRange[i]
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


function Simulator(boxSizes::Vector{Int64}, inputGridVectors::Matrix{Float64}, parameters::Parameters)
    log_info("🚀 Initializing the simulator...")
    box = CreateBoxByPrimaryVectors(parameters.primaryVectors, boxSizes)
    ranges = parameters.latticeRanges 
    atomNumber = (ranges[1,2] - ranges[1,1]) * (ranges[2,2] - ranges[2,1]) * (ranges[3,2] - ranges[3,1]) * length(parameters.basisTypes)
    log_info("🔬 Number of atoms in the box: $(atomNumber)  (= $(ranges[1,2] - ranges[1,1]) × $(ranges[2,2] - ranges[2,1]) × $(ranges[3,2] - ranges[3,1]) × $(length(parameters.basisTypes)))\n")
    simulator = Simulator(box, inputGridVectors, parameters)
    if !IS_DYNAMIC_LOAD
        _initSimulatorAtoms!(simulator, parameters)
    end
    log_success("🎉 Simulator initialized!\n")
    return simulator    
end 

function Simulator(boxVectors::Matrix{Float64}, inputGridVectors::Matrix{Float64}, parameters::Parameters)
    # this function should be improved by warining is orthogonal and type error when desird for boxsize but in float
    log_info("🚀 Initializing the simulator...")
    box = Box(boxVectors)
    ranges = parameters.latticeRanges 
    atomNumber = (ranges[1,2] - ranges[1,1]) * (ranges[2,2] - ranges[2,1]) * (ranges[3,2] - ranges[3,1]) * length(parameters.basisTypes)
    log_info("🔬 Number of atoms in the box: $(atomNumber)  (= $(ranges[1,2] - ranges[1,1]) × $(ranges[2,2] - ranges[2,1]) × $(ranges[3,2] - ranges[3,1]) × $(length(parameters.basisTypes)))\n")
    simulator = Simulator(box, inputGridVectors, parameters)
    if !IS_DYNAMIC_LOAD
        _initSimulatorAtoms!(simulator, parameters)
    end
    log_success("🎉 Simulator initialized!\n")
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
    log_success("🧪 $(simulator.numberOfAtoms) atoms created.\n")
    InitLatticePointEnvronment(simulator)
    for cell in simulator.grid.cells
        cell.atomicDensity = length(cell.atoms) / simulator.grid.cellVolume
    end 
    for atom in simulator.atoms
        Pertubation!(atom, simulator)
    end
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


function DisplaceAtom!(atom::Atom, newPosition::Vector{Float64}, simulator::Simulator)
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
    #@show newPosition, atom.coordinate
    SetCoordinate!(atom, newPosition)
    cellIndex = WhichCell(atom.coordinate, simulator.grid)

    idx = LinearIndex(atom.cellIndex, simulator.grid.sizes)
    if cellIndex != atom.cellIndex
        ChangeCell!(atom, cellIndex, simulator)
    end
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


function ComputeP!(atom_p::Atom, atom_t::Atom, crossFlag::NTuple{3, Int8}, box::Box, pMax::Float64)
    dv = VectorDifference(atom_p.coordinate, atom_t.coordinate, crossFlag, box)
    #@show dv
    #for d in dv
    #    if abs(d) > pMax
    #        return Inf
    #    end
    #end
    t = dot(dv, atom_p.velocityDirection)
    atom_t.pL = t
    atom_t.pPoint = atom_p.coordinate + t * atom_p.velocityDirection
    atom_t.pVector = atom_t.pPoint - atom_t.coordinate
    p = norm(atom_t.pVector)
    atom_t.pValue = p
    # need to check periodic condition
    #if atom_t.index == 3055
    #    @show atom_p.coordinate, atom_t.coordinate, p, atom_t.pValue, dv, atom_t.pPoint, atom_t.pVector
    #end
    return p
end




function SimultaneousCriteria(atom::Atom, neighborAtom::Atom, addedTarget::Atom, simulator::Simulator)
    flagP = abs(neighborAtom.pValue - addedTarget.pValue) <= simulator.parameters.pMax
    flagQ = abs(neighborAtom.pL - addedTarget.pL) <= simulator.constantsByType.qMax[[neighborAtom.type, addedTarget.type]]
    return flagP && flagQ
end


function ChangeCell!(atom::Atom, nextCellIndex::Tuple{Int64, Int64, Int64}, simulator::Simulator)
    originalCell = GetCell(simulator.grid, atom.cellIndex)
    delete!(originalCell, atom, simulator)
    nextCell = GetCell(simulator.grid, nextCellIndex)
    push!(nextCell, atom, simulator)
end


function SetVelocityDirection!(atom::Atom, velocity::Vector{Float64})
    n = norm(velocity)
    if isnan(n) || n == Inf
        atom.velocityDirection = [0.0, 0.0, 0.0]
    else
        atom.velocityDirection = velocity / n
    end
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
        index = neighborCellInfo.index
        cross = neighborCellInfo.cross
        neighborCell = GetCell(grid, index)
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
    if simulator.parameters.temperature > 0.0
        for d in 1:3
            atom.coordinate[d] += GaussianDeltaX(simulator.constantsByType.sigma[atom.type])
        end
    end
end

function SetCoordinate!(atom::Atom, coordinate::Vector{Float64})
    atom.coordinate .= coordinate
end 


function TemperatureToSigma(T::Float64, θ_D::Float64, m_rel::Float64; atol=1e-10, rtol=1e-8)
    if T == 0.0
        log_info("🌡️ Temperature is 0.0 K")
        return 0
    end
    ħ   = 1.054_571_817e-34      # J·s
    kB  = 1.380_649_000e-23      # J/K
    amu = 1.660_539_066_60e-27   # kg

    M = m_rel * amu
    y_max = θ_D / T      

    # 积分 ∫ x/(e^x-1) dx
    integrand(x) = x / (exp(x) - 1)
    I, _ = quadgk(integrand, 0.0, y_max; atol, rtol)

    σ2 = 3 * ħ^2 / (M * kB * θ_D) * (0.25 + (T/θ_D)^2 * I)
    σ  = sqrt(σ2) * 1e10         # m → Å

    return σ
end

