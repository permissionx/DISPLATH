using StaticArrays

function InitCellStd!(simulator::Simulator, primaryCellNumbersInCell::Vector{Int64})
    nP = primaryCellNumbersInCell
    parameters = simulator.parameters
    basisTypes = parameters.basisTypes
    basis = parameters.basis
    primaryVectors = parameters.primaryVectors
    cellsStd = simulator.cellsStd
    for _ in 1:2
        indexInCell = 0
        cellStd = CellStd()
        for X in 0:nP[1]-1
            for Y in 0:nP[2]-1
                for Z in 0:nP[3]-1
                    for i in 1:length(basisTypes)
                        indexInCell += 1
                        x = (X + basis[i,1]) * primaryVectors[1,1]
                        y = (Y + basis[i,2]) * primaryVectors[2,2]
                        z = (Z + basis[i,3]) * primaryVectors[3,3]
                        atom = Atom(basisTypes[i], [x, y, z], parameters)
                        atom.index = 0 
                        atom.indexInCell = indexInCell
                        push!(cellStd.atoms, atom)
                    end
                end
            end
        end
        push!(cellsStd, cellStd)
    end
end


function ComputeLatticeAtoms_Orthogonal!(cell::Cell, simulator::Simulator)
    cellsStd = simulator.cellsStd
    vacancyIDs = [v.indexInCell for v in cell.vacancies]
    nLatticeAtoms = 0
    latticeRanges = simulator.parameters.latticeRanges
    index = cell.index
    isEmpty = false
    for d in 1:3
        if index[d] < latticeRanges[d,1] / simulator.parameters.primaryCellNumbersInCell[d] || index[d] > latticeRanges[d,2] / simulator.parameters.primaryCellNumbersInCell[d]
            isEmpty = true
        end
    end
    for i in 1:length(cellsStd[1].atoms)
        if isEmpty || i in vacancyIDs
            cellsStd[1].atoms[i].isAlive = false
        else
            nLatticeAtoms += 1
            for d in 1:3
                cellsStd[1].atoms[i].coordinate[d] = cellsStd[2].atoms[i].coordinate[d] + cell.ranges[d,1]
            end
            cellsStd[1].atoms[i].isAlive = true
            Pertubation!(cellsStd[1].atoms[i], simulator)
        end
    end
    cell.atomicDensity = (length(cell.atoms) + nLatticeAtoms) / simulator.grid.cellVolume
end




function LoadCellAtoms!(cell::Cell, simulator::Simulator)
    ComputeLatticeAtoms_Orthogonal!(cell, simulator)
end

function SyncLatticeAtomBuffer(atom::Atom, cell::Cell, simulator::Simulator)
    # from stdCell to buff
    # atom is from stdCell
    count = simulator.latticeTargetsBuffer.count
    buffAtom = simulator.latticeTargetsBuffer.atoms[count]
    buffAtom.type = atom.type
    buffAtom.coordinate[:] = atom.coordinate[:]
    buffAtom.indexInCell = atom.indexInCell
    buffAtom.cellIndex = cell.index  # need to check if correct
    buffAtom.pL = atom.pL
    buffAtom.pPoint = SVector{3, Float64}(atom.pPoint[1], atom.pPoint[2], atom.pPoint[3])
    buffAtom.pVector = SVector{3, Float64}(atom.pVector[1], atom.pVector[2], atom.pVector[3])
    buffAtom.pValue = atom.pValue
    simulator.latticeTargetsBuffer.count += 1
    return buffAtom
end




function GetTargetsFromNeighbor_dynamicLoad(atom::Atom, cell::Cell, filterIndexes::Vector{Int64}, filterLatticeIndexes::Vector{Tuple{Int64,Int64,Int64,Int64}}, simulator::Simulator)
    grid = simulator.grid
    box = simulator.box
    targets = Vector{Atom}()
    pMax = simulator.parameters.pMax
    nthreads = Threads.nthreads()
    neighborCellsInfo = cell.neighborCellsInfo
    AlreadyLoadedFlags = [true for _ in 1:27]
    infiniteFlag_tls = [true for _ in 1:27]
    threadCandidates = simulator.workBuffers.threadCandidates
    latticeRanges = simulator.parameters.latticeRanges  # need to exclude no atom cell
    for tc in threadCandidates
        empty!(tc)
    end
    for n in 1:length(neighborCellsInfo)
        neighborCellInfo = neighborCellsInfo[n]
        index = neighborCellInfo.index
        GetCell(grid, index)  # preload the cell to avoid race condition
    end
    for n in 1:length(neighborCellsInfo)
        neighborCellInfo = neighborCellsInfo[n]
        cross = neighborCellInfo.cross
        nonPeriodicFlag = false 
        for d in 1:3
            if cross[d] != 0 && !simulator.parameters.periodic[d] 
                nonPeriodicFlag = true
                break
            end
        end
        if nonPeriodicFlag
            continue 
        end 
        buf = threadCandidates[Threads.threadid()]
        index = neighborCellInfo.index
        neighborCell = GetCell(grid, index)
        if neighborCell.isExplored 
            continue
        end
        LoadCellAtoms!(neighborCell, simulator)
        neighborCell.isExplored = true
        infiniteFlag_tls[n] = false
        na = length(neighborCell.atoms) 
        
        #for neighborAtom in [neighborCell.atoms; simulator.cellsStd[1].atoms]
        for n in 1:na + length(simulator.cellsStd[1].atoms)
            if n <= na
                neighborAtom = neighborCell.atoms[n]
            else
                neighborAtom = simulator.cellsStd[1].atoms[n-na]
            end
            latticeIndex = (neighborAtom.cellIndex[1], neighborAtom.cellIndex[2], neighborAtom.cellIndex[3], neighborAtom.indexInCell)
            if (n <= na && (neighborAtom.index == atom.index || neighborAtom.index in filterIndexes))  || (n > na && (!neighborAtom.isAlive || latticeIndex in filterLatticeIndexes))    
                continue
            end
            if ComputeVDistance(atom, neighborAtom, neighborCellInfo.cross, box) > 0 
                p = ComputeP!(atom, neighborAtom, neighborCellInfo.cross, box)
                if p >= pMax
                    continue
                end
                if n > na
                    neighborAtom = SyncLatticeAtomBuffer(neighborAtom, neighborCell, simulator)
                end
                push!(buf, neighborAtom)
            end
        end
    end

    candidateTargets = simulator.workBuffers.candidateTargets
    empty!(candidateTargets)
    for tc in threadCandidates
        append!(candidateTargets, tc)
    end
    infiniteFlag = reduce(&, infiniteFlag_tls)
    for neighborCellInfo in neighborCellsInfo
        idx = neighborCellInfo.index
        cell = GetCell(simulator.grid, idx)
        push!(simulator.exploredCells, cell)
    end

    if isempty(candidateTargets)
        return (targets, infiniteFlag)
    end
    _, minIdx = findmin(t -> t.pL, candidateTargets)
    nearestTarget = candidateTargets[minIdx]   
    push!(targets, nearestTarget)


    for candidateTarget in candidateTargets
        if candidateTarget.index == nearestTarget.index
            continue
        end
        if SimultaneousCriteria(candidateTarget, nearestTarget, simulator)
            push!(targets, candidateTarget)
        end
    end    
    return (targets, infiniteFlag)
end

function SimultaneousCriteria(candidateTarget::Atom, nearestTarget::Atom, simulator::Simulator)
    deltaPL = candidateTarget.pL - nearestTarget.pL
    if deltaPL > simulator.constantsByType.qMax[[candidateTarget.type, nearestTarget.type]]
        return false
    elseif nearestTarget.pValue * nearestTarget.pValue + deltaPL * deltaPL > simulator.parameters.pMax_squared 
        return false
    elseif candidateTarget.pValue * candidateTarget.pValue + deltaPL * deltaPL > simulator.parameters.pMax_squared 
        return false
    end
    return true
end

function Collision_dynamicLoad!(atom_p::Atom, atoms_t::Vector{Atom}, simulator::Simulator)
    N_t = length(atoms_t)
    grid = simulator.grid
    buffers = simulator.workBuffers.collisionParames
    EnsureCollisionCapacity!(buffers, N_t)
    tanφList = @view buffers.tanφList[1:N_t]
    tanψList = @view buffers.tanψList[1:N_t]
    E_tList = @view buffers.E_tList[1:N_t]
    x_pList = @view buffers.x_pList[1:N_t]
    x_tList = @view buffers.x_tList[1:N_t]
    Q_locList = @view buffers.Q_locList[1:N_t]
    atom_t = atoms_t[1]
    pL = atom_t.pL   
    pPoint = atom_t.pPoint
    pL -= atom_p.emptyPath
    N = simulator.uniformDensity
    Q_nl_v = Q_nl(atom_p.energy, atom_p.mass, atom_t.mass, atom_p.type, atom_t.type,
                         pL, N, simulator.constantsByType)
    atom_p.energy -= Q_nl_v
    if atom_p.energy < 0.1 && atom_p.energy + Q_nl_v >= 0.1
        atom_p.energy = 0.11
    end
    momentum = @SVector [0.0, 0.0, 0.0] 
    for (i, atom_t) in enumerate(atoms_t)
        p = atom_t.pValue
        #N = simulator.uniformDensity 
        tanφList[i], tanψList[i], E_tList[i], x_pList[i], x_tList[i], Q_locList[i] = CollisionParams(
            atom_p.energy, atom_p.mass, atom_t.mass, atom_p.type, atom_t.type, p, simulator.constantsByType,
            simulator.θFunctions[[atom_p.type, atom_t.type]], simulator.τFunctions[[atom_p.type, atom_t.type]])   
        if atom_t.pValue != 0
            velocityDirectionTmp = -atom_t.pVector / atom_t.pValue * tanψList[i] + atom_p.velocityDirection
        else
            velocityDirectionTmp = atom_p.velocityDirection
        end   
        SetVelocityDirection!(atom_t, velocityDirectionTmp)
        momentum += sqrt(2 * atom_t.mass * E_tList[i]) * atom_t.velocityDirection
    end
    pMomentum = sqrt(2 * atom_p.mass * atom_p.energy) * atom_p.velocityDirection - momentum
    pVelocity = pMomentum  / atom_p.mass
    SetVelocityDirection!(atom_p, pVelocity)
    pEnergy =  sum(pMomentum .* pMomentum) / 2 / atom_p.mass
    sumE_t = sum(E_tList)
    sumQ_loc = sum(Q_locList) 
    ENeed = atom_p.energy - sumQ_loc # - (N_t - 1) * Q_nl_v
    λ = ENeed / (pEnergy + sumE_t)
    DisplaceAtom!(atom_p, pPoint, simulator)
    SetEnergy!(atom_p, pEnergy * λ)
    #if atom_p.type == 2
    #    @record "log/$(simulator.nCascade).csv" "$(pEnergy * λ),$(minimum([a.pValue for a in atoms_t])),$(pL),$(N_t),$(atom_p.coordinate[1]),$(atom_p.coordinate[2]),$(atom_p.coordinate[3]),$(atom_p.velocityDirection[1]),$(atom_p.velocityDirection[2]),$(atom_p.velocityDirection[3])" "e,p,pL,N_t,x,y,z,vx,vy,vz" 
    #end
    E_tList *= λ
    for (i, atom_t) in enumerate(atoms_t)
        if E_tList[i] > GetDTE(atom_t, simulator) && E_tList[i] - GetBDE(atom_t, simulator) > 0.1
            SetEnergy!(atom_t, E_tList[i] - GetBDE(atom_t, simulator))
        else
            SetEnergy!(atom_t, 0.0)
        end
    end
end 


function DumpInCascade_dynamicLoad(simulator::Simulator)
    if simulator.parameters.isDumpInCascade
        if simulator.parameters.debugMode == false
            @dump "Cascade_$(simulator.nCascade).dump" [simulator.atoms; simulator.vacancies] ["vx", "vy", "vz", "e"]
        else
            cells = values(simulator.grid.cells)
            b = [atom for cell in cells for atom in cell.atoms]
            @dump "Cascade_$(simulator.nCascade).dump" b ["vx", "vy", "vz", "e"]
        end
    end
end


function Cascade_dynamicLoad!(atom_p::Atom, simulator::Simulator)
    pAtoms = Vector{Atom}([atom_p])
    pAtomsIndex = [a.index for a in pAtoms]
    parameters = simulator.parameters
    simulator.nCollisionEvent = 0
    simulator.nCascade += 1
    DumpInCascade_dynamicLoad(simulator)
    while true
        simulator.nCollisionEvent += 1
        targetsList = Vector{Vector{Atom}}()
        deleteIndexes = Int64[]
        othersTargetIndexes = Int64[]
        othersLatticeTargetIndexes = Vector{Tuple{Int64, Int64, Int64, Int64}}()
        simulator.latticeTargetsBuffer.count = 1
        for (na, pAtom) in enumerate(pAtoms)
            targets, isAlive = ShotTarget_dynamicLoad(pAtom, [pAtomsIndex; pAtom.lastTargets; othersTargetIndexes],othersLatticeTargetIndexes, simulator)
            if !isAlive
                empty!(pAtom.lastTargets)
                delete_dynamicLoad!(simulator, pAtom)
                push!(deleteIndexes, na)
                continue
            end
            push!(targetsList, targets)
            for t in targets
                if t.isNewlyLoaded
                    push!(othersLatticeTargetIndexes, (t.cellIndex[1], t.cellIndex[2], t.cellIndex[3], t.indexInCell))
                else
                    push!(othersTargetIndexes, t.index)
                end
            end
        end
        deleteat!(pAtoms, deleteIndexes)
        pAtomsIndex = [a.index for a in pAtoms]
        nextPAtoms = Vector{Atom}()
        for (pAtom, targets) in zip(pAtoms, targetsList)
            if length(targets) > 0
                pAtom.lastTargets = [t.index for t in targets]
                Collision_dynamicLoad!(pAtom, targets, simulator)
                for target in targets
                    if target.energy > 0.0   
                        targetEntity = LeaveLatticePoint_dynamicLoad!(target, simulator)
                        DisplaceAtom!(targetEntity, targetEntity.coordinate, simulator)
                        push!(nextPAtoms, targetEntity)
                        targetEntity.lastTargets = [pAtom.index]
                    end
                end
                if pAtom.energy > parameters.stopEnergy 
                    push!(nextPAtoms, pAtom)
                else
                    pAtom.lastTargets = Vector{Int64}()
                    Stop_dynamicLoad!(pAtom, simulator)
                end
            else
                push!(nextPAtoms, pAtom)
            end
        end
        DumpInCascade_dynamicLoad(simulator)
        if length(nextPAtoms) > 0
            pAtoms = nextPAtoms
            sort!(pAtoms, by = a -> a.energy, rev = true)
            pAtomsIndex = [a.index for a in pAtoms]
        else
            break
        end
    end
end


function delete_dynamicLoad!(simulator::Simulator, atom::Atom; isDeleteVacancy::Bool = false)
    cell = GetCell(simulator.grid, atom.cellIndex)
    if !isDeleteVacancy
        deleteat!(cell.atoms, findfirst(a -> a.index == atom.index, cell.atoms))
        simulator.numberOfAtoms -= 1
    else
        deleteat!(cell.vacancies, findfirst(v -> v.index == atom.index, cell.vacancies))
        simulator.numberOfVacancies -= 1
    end
    atom.isAlive = false 
end

function Stop_dynamicLoad!(atom::Atom, simulator::Simulator)
    grid = simulator.grid
    cell = GetCell(grid, atom.cellIndex)
    nearestVacancyDistance_squared = Inf
    isExist = false
    nearestVacancy = nothing  
    if ! cell.isPushedNeighbor
        SetCellNeighborInfo!(cell, grid)
        cell.isPushedNeighbor = true
    end
    if simulator.parameters.vacancyRecoverDistance_squared == 0.0
        return
    end
    for neighborCellInfo in cell.neighborCellsInfo
        index = neighborCellInfo.index
        cross = neighborCellInfo.cross
        neighborCell = GetCell(simulator.grid, index)
        for vacancy in neighborCell.vacancies
            dr2 = ComputeDistance_squared(atom.coordinate, vacancy.coordinate, cross, simulator.box)
            if dr2 < simulator.parameters.vacancyRecoverDistance_squared && dr2 < nearestVacancyDistance_squared
                nearestVacancyDistance_squared = dr2
                nearestVacancy = vacancy  # store the nearest vacancy
                isExist = true
            end
        end
    end
    if isExist && nearestVacancy !== nothing
        if atom.type == nearestVacancy.type - length(keys(simulator.parameters.typeDict)) 
            delete_dynamicLoad!(simulator, atom)
            delete_dynamicLoad!(simulator, nearestVacancy, isDeleteVacancy = true)
        else
            SetCoordinate!(atom, nearestVacancy.coordinate)
            Pertubation!(atom, simulator)
            ChangeCell!(atom, nearestVacancy.cellIndex, simulator)
        end
    end
end


function LeaveLatticePoint_dynamicLoad!(buffAtom::Atom, simulator::Simulator; isUpdateEnv::Bool = true)
    if buffAtom.isNewlyLoaded
        atom = CopyAtom(buffAtom, simulator)
        cell = GetCell(simulator.grid, atom.cellIndex)
        vacancy = CreateVacancy(atom, simulator)
        push!(cell.vacancies, vacancy)
        push!(simulator.vacancies, vacancy)
        simulator.numberOfVacancies += 1
        vacancy.index = simulator.maxVacancyID
        simulator.maxVacancyID += 1


        push!(cell.atoms, atom)
        push!(simulator.atoms, atom)
        simulator.maxAtomID += 1
        atom.index = simulator.maxAtomID
        simulator.numberOfAtoms += 1
        return atom
    else
        return buffAtom
    end
end

function CopyAtom(atom::Atom, simulator::Simulator)
    coord = if atom.coordinate isa SVector
        [atom.coordinate[1], atom.coordinate[2], atom.coordinate[3]]
    else
        atom.coordinate[:]
    end
    newAtom = Atom(atom.type, coord, simulator.parameters)
    newAtom.cellIndex = atom.cellIndex
    newAtom.pL = atom.pL
    newAtom.pPoint = SVector{3, Float64}(atom.pPoint[1], atom.pPoint[2], atom.pPoint[3])
    newAtom.pVector = SVector{3, Float64}(atom.pVector[1], atom.pVector[2], atom.pVector[3])
    newAtom.pValue = atom.pValue
    return newAtom
end

function CreateVacancy(atom::Atom, simulator::Simulator)
    coord = if atom.coordinate isa SVector
        [atom.coordinate[1], atom.coordinate[2], atom.coordinate[3]]
    else
        atom.coordinate[:]
    end
    vacancy = Atom(atom.type, coord, simulator.parameters)
    vacancy.cellIndex = atom.cellIndex
    vacancy.type += length(keys(simulator.parameters.typeDict))
    vacancy.indexInCell = atom.indexInCell
    return vacancy
end


                                                   


function ShotTarget_dynamicLoad(atom::Atom, filterIndexes::Vector{Int64}, filterLatticeIndexes::Vector{Tuple{Int64,Int64,Int64,Int64}}, simulator::Simulator)
    grid = simulator.grid
    periodic = simulator.parameters.periodic    
    cell = GetCell(grid, atom.cellIndex)
    atom.emptyPath = 0.0
    while true
        if ! cell.isPushedNeighbor
            SetCellNeighborInfo!(cell, grid)
            cell.isPushedNeighbor = true
        end
        targets, isInfinity = GetTargetsFromNeighbor_dynamicLoad(atom, cell, filterIndexes, filterLatticeIndexes, simulator)
        if length(targets) > 0
            for cell in simulator.exploredCells
                cell.isExplored = false
            end
            empty!(simulator.exploredCells)
            return targets, true
        else
            dimension, direction, t = AtomOutFaceDimension(atom, cell)
            atom.emptyPath = t
            neighborIndex = MVector{3,Int8}(0, 0, 0)  
            neighborIndex[dimension] = direction == 1 ? Int8(-1) : Int8(1)
            neighborIndex .+= 2
            neighborInfo = cell.neighborCellsInfo[neighborIndex[1], neighborIndex[2], neighborIndex[3]]
            crossFlag = neighborInfo.cross
            if crossFlag[dimension] != 0 && periodic[dimension]
                atom.coordinate[dimension] -= crossFlag[dimension] * simulator.box.vectors[dimension, dimension]
            end
            if (neighborInfo.cross[dimension] != 0 && !periodic[dimension]) || isInfinity
                for cell in simulator.exploredCells
                    cell.isExplored = false
                end
                atom.emptyPath = 0.0
                empty!(simulator.exploredCells)
                return Vector{Atom}(), false # means find nothing  
            end 
            index = neighborInfo.index
            cell = GetCell(grid, index)
        end
    end
end


function Dump_dynamicLoad(simulator::Simulator, fileName::String, step::Int64, type::String="a", isDebug::Bool=false)
    if !simulator.parameters.isOrthogonal        
        error("The box is not orthogonal, please use the orthogonal box.")
    end
    open(fileName, type) do file
        write(file, "ITEM: TIMESTEP\n")
        write(file, string(step), "\n")
        write(file, "ITEM: NUMBER OF ATOMS\n")
        write(file, string(simulator.numberOfAtoms + simulator.numberOfVacancies), "\n")
        write(file, "ITEM: BOX BOUNDS ")
        for d in 1:3
            if simulator.parameters.periodic[d] 
                write(file, "pp ")
            else
                write(file, "ff ")
            end
        end
        write(file, "\n")
        for d in 1:3
            write(file, "0 $(simulator.box.vectors[d,d])\n")
        end
        if isDebug
            write(file, "ITEM: ATOMS id type x y z vx vy vz energy cx cy cz dte\n")
        else
            write(file, "ITEM: ATOMS id type x y z e\n")
        end
        for atom in simulator.atoms
            if atom.isAlive
                if isDebug
                    write(file, "$(atom.index) $(atom.type) \
                    $(atom.coordinate[1]) $(atom.coordinate[2]) $(atom.coordinate[3]) \
                    $(atom.velocityDirection[1]*sqrt(2*atom.mass*atom.energy)) $(atom.velocityDirection[2]*sqrt(2*atom.mass*atom.energy)) $(atom.velocityDirection[3]*sqrt(2*atom.mass*atom.energy)) \
                    $(atom.energy) \
                    $(atom.cellIndex[1]) $(atom.cellIndex[2]) $(atom.cellIndex[3]) \
                    $(GetDTE(atom, simulator))\n")
                else
                    write(file, "$(atom.index) $(atom.type) \
                    $(atom.coordinate[1]) $(atom.coordinate[2]) $(atom.coordinate[3]) $(atom.energy)\n")
                end
            end
        end 
        for atom in simulator.vacancies
            if atom.isAlive
                if isDebug
                    write(file, "$(atom.index+100000) $(atom.type) \
                    $(atom.coordinate[1]) $(atom.coordinate[2]) $(atom.coordinate[3]) \
                    $(atom.velocityDirection[1]*sqrt(2*atom.mass*atom.energy)) $(atom.velocityDirection[2]*sqrt(2*atom.mass*atom.energy)) $(atom.velocityDirection[3]*sqrt(2*atom.mass*atom.energy)) \
                    $(atom.energy) \
                    $(atom.cellIndex[1]) $(atom.cellIndex[2]) $(atom.cellIndex[3]) \
                    $(GetDTE(atom, simulator))\n")
                else
                    write(file, "$(atom.index+100000) $(atom.type) \
                    $(atom.coordinate[1]) $(atom.coordinate[2]) $(atom.coordinate[3]) $(atom.energy)\n")
                end
            end
        end
    end
end



function Restore_dynamicLoad!(simulator::Simulator)
    parameters = simulator.parameters
    for atom in [simulator.atoms; simulator.vacancies]
        if atom.isAlive
            cellIndex = atom.cellIndex
            cell = GetCell(simulator.grid, cellIndex)
            empty!(cell.atoms)
        end
    end
    empty!(simulator.atoms)
    empty!(simulator.vacancies)
    simulator.maxAtomID = 0
    simulator.maxVacancyID = 1E6 
    simulator.numberOfAtoms = 0
    simulator.numberOfVacancies = 0
end

