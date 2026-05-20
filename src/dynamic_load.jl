using StaticArrays

function InitCellStd!(simulator::Simulator, PN::Vector{Int64})
    parameters = simulator.parameters
    basisTypes = parameters.basisTypes
    basis = parameters.basis
    primaryVectors = parameters.primaryVectors
    cellsStd = simulator.cellsStd
    indexInCell = 0
    cellStd = simulator.cellStd
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
                    push!(cellStd.atoms, atprimaryVectorsom)
                end
            end
        end
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
        atom = cellsStd[1].atoms[i]
        if isEmpty || i in vacancyIDs
            atom.isAlive = false
        else
            nLatticeAtoms += 1
            for d in 1:3
                atom.coordinate[d] = cellsStd[2].atoms[i].coordinate[d] + cell.ranges[d,1]
            end
            atom.isAlive = true
            atom.cellIndex = cell.index
            Pertubation!(atom, simulator)
        end
    end
    cell.atomicDensity = (length(cell.atoms) + nLatticeAtoms) / simulator.grid.cellVolume
end


function LoadCellAtoms!(cell::Cell, simulator::Simulator)
        if simulator.parameters.isPrimaryVectorOrthogonal
            ComputeLatticeAtoms_Orthogonal!(cell, simulator)
        else
            error("Lattice must be orthogonal!")
        end
        cell.atomicDensity = length(cell.latticeAtoms) / simulator.grid.cellVolume
end

function GetTargetsFromNeighbor_dynamicLoad(atom::Atom, cell::Cell, filterIndexes::Vector{Int64}, simulator::Simulator)
    grid = simulator.grid
    box = simulator.box
    targets = Vector{Atom}()
    pMax = simulator.parameters.pMax
    nthreads = Threads.nthreads()
    neighborCellsInfo = cell.neighborCellsInfo
    threadCandidates = simulator.workBuffers.threadCandidates
    for tc in threadCandidates
        empty!(tc)
    end
    if nthreads >= 1 
        for n in 1:length(neighborCellsInfo)
            neighborCellInfo = neighborCellsInfo[n]
            index = neighborCellInfo.index
            GetCell(grid, index)  # preload the cell to avoid race condition
        end
    end
    @threads :static for n in 1:length(neighborCellsInfo)
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
        LoadCellAtoms!(neighborCell, simulator)
        na = length(neighborCell.atoms) 
        for n in 1:na+ length(neighborCell.latticeAtoms)
        #for neighborAtom in [neighborCell.atoms; neighborCell.latticeAtoms]
            if n <= na
                neighborAtom = neighborCell.atoms[n]
            else
                neighborAtom = neighborCell.latticeAtoms[n - na]
            end
            if neighborAtom.index == atom.index || neighborAtom.index in filterIndexes    
                continue
            end
            if ComputeVDistance(atom, neighborAtom, neighborCellInfo.cross, box) > 0 
                p = ComputeP!(atom, neighborAtom, neighborCellInfo.cross, box)
                if p >= pMax
                    continue
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


    if isempty(candidateTargets)
        return targets
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
    return targets
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
    #if simulator.nCollisionEvent < 3
    #    @show atom_p.emptyPath
    #    @show pL
    #    @show N
    #    @show simulator.uniformDensity
    #elseif simulator.nCollisionEvent == 3
    #    exit()
    #end
    Q_nl_v = Q_nl(atom_p.energy, atom_p.mass, atom_t.mass, atom_p.type, atom_t.type,
                         pL, N, simulator.constantsByType)
    atom_p.energy -= Q_nl_v
    #if atom_p.type == 2
    #global Q_loss += Q_nl_v  # debug 
    #end
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
            a = [atom for cell in cells for atom in cell.latticeAtoms]
            b = [atom for cell in cells for atom in cell.atoms]
            @dump "Cascade_$(simulator.nCascade).dump" [a; b] ["vx", "vy", "vz", "e"]
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
        for (na, pAtom) in enumerate(pAtoms)
            targets, isAlive = ShotTarget_dynamicLoad(pAtom, [pAtomsIndex; pAtom.lastTargets; othersTargetIndexes], simulator)
            if !isAlive
                empty!(pAtom.lastTargets)
                delete_dynamicLoad!(simulator, pAtom)
                push!(deleteIndexes, na)
                continue
            end
            push!(targetsList, targets)
            append!(othersTargetIndexes, [t.index for t in targets])
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
                        LeaveLatticePoint_dynamicLoad!(target, simulator)
                        DisplaceAtom!(target, target.coordinate, simulator)
                        push!(nextPAtoms, target)
                        target.lastTargets = [pAtom.index]
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
    if simulator.nCascade % parameters.nCascadeEveryLoad == 0
        rss = parse(Int, read(`ps -o rss= -p $(getpid())`, String)) 
        if rss > simulator.parameters.maxRSS
            CleanUpLatticeAtoms(simulator)
        end
    end
end



function CleanUpLatticeAtoms(simulator::Simulator)
    empty!(simulator.grid.cells)
    GC.gc()
    for atom in simulator.atoms
        if atom.isAlive 
            cell = GetCell(simulator.grid, atom.cellIndex)
            push!(cell.atoms, atom)
        end
    end
    for vacancy in simulator.vacancies
        if vacancy.isAlive
            cell = GetCell(simulator.grid, vacancy.cellIndex)
            push!(cell.vacancies, vacancy)
        end
    end
    simulator.minLatticeAtomID = 0
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


function LeaveLatticePoint_dynamicLoad!(atom::Atom, simulator::Simulator; isUpdateEnv::Bool = true)
    if atom.isNewlyLoaded
        cell = GetCell(simulator.grid, atom.cellIndex)
        atom.isNewlyLoaded = false
        vacancy = CreateVacancy(atom, simulator)
        push!(cell.vacancies, vacancy)
        push!(simulator.vacancies, vacancy)
        simulator.numberOfVacancies += 1
        vacancy.index = simulator.maxVacancyID
        simulator.maxVacancyID += 1

        deleteat!(cell.latticeAtoms, findfirst(a -> a.index == atom.index, cell.latticeAtoms))
        push!(cell.atoms, atom)
        push!(simulator.atoms, atom)
        simulator.maxAtomID += 1
        atom.index = simulator.maxAtomID
        simulator.numberOfAtoms += 1
    end
end

function CopyAtom(atom::Atom, simulator::Simulator)
    newAtom = Atom(atom.type, atom.coordinate, simulator.parameters)
    newAtom.cellIndex = atom.cellIndex
    return newAtom
end

function CreateVacancy(atom::Atom, simulator::Simulator)
    coord = if atom.latticeCoordinate isa SVector
        [atom.latticeCoordinate[1], atom.latticeCoordinate[2], atom.latticeCoordinate[3]]
    else
        atom.latticeCoordinate[:]
    end
    vacancy = Atom(atom.type, coord, simulator.parameters)
    vacancy.cellIndex = atom.cellIndex
    vacancy.type += length(keys(simulator.parameters.typeDict))
    vacancy.indexInCell = atom.indexInCell
    return vacancy
end


                                                   


function ShotTarget_dynamicLoad(atom::Atom, filterIndexes::Vector{Int64}, simulator::Simulator)
    grid = simulator.grid
    periodic = simulator.parameters.periodic    
    cell = GetCell(grid, atom.cellIndex)
    atom.emptyPath = 0.0
    while true
        if ! cell.isPushedNeighbor
            SetCellNeighborInfo!(cell, grid)
            cell.isPushedNeighbor = true
        end
        targets = GetTargetsFromNeighbor_dynamicLoad(atom, cell, filterIndexes, simulator)
        if length(targets) > 0
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
            if (neighborInfo.cross[dimension] != 0 && !periodic[dimension]) || t >= simulator.parameters.infiniteLength
                atom.emptyPath = 0.0
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

