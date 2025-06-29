function ComputeLatticeAtoms_Orthogonal!(cell::Cell, simulator::Simulator)
    parameters = simulator.parameters
    primaryVectors = parameters.primaryVectors
    latticeRanges = parameters.latticeRanges
    basisTypes = parameters.basisTypes
    basis = parameters.basis
    if !cell.isSavedAtomRange
        # Optimized calculation for orthogonal case
        vertexMatrix = cell.vertexMatrix
        a1, a2, a3 = primaryVectors[1,1], primaryVectors[2,2], primaryVectors[3,3]
        # Direct calculation without matrix multiplication
        nfrac1 = vertexMatrix[:, 1] ./ a1
        nfrac2 = vertexMatrix[:, 2] ./ a2  
        nfrac3 = vertexMatrix[:, 3] ./ a3
        nmin = [floor(Int, minimum(nfrac1) - 1), 
                floor(Int, minimum(nfrac2) - 1), 
                floor(Int, minimum(nfrac3) - 1)]
        nmax = [ceil(Int, maximum(nfrac1) + 1),
                ceil(Int, maximum(nfrac2) + 1),
                ceil(Int, maximum(nfrac3) + 1)]
        n1 = max(nmin[1],latticeRanges[1,1]):min(nmax[1],latticeRanges[1,2])
        n2 = max(nmin[2],latticeRanges[2,1]):min(nmax[2],latticeRanges[2,2])
        n3 = max(nmin[3],latticeRanges[3,1]):min(nmax[3],latticeRanges[3,2])
        cell.atomRange = [n1, n2, n3]
        cell.isSavedAtomRange = true
    else
        n1, n2, n3 = cell.atomRange
    end
    # Get cell ranges for coordinate filtering
    ranges = cell.ranges
    x_min, x_max = ranges[1, 1], ranges[1, 2]
    y_min, y_max = ranges[2, 1], ranges[2, 2] 
    z_min, z_max = ranges[3, 1], ranges[3, 2]
    atoms = Vector{Atom}()
    # Optimized generation for orthogonal case
    a1, a2, a3 = primaryVectors[1,1], primaryVectors[2,2], primaryVectors[3,3]
    # Pre-compute which lattice ranges actually intersect with cell
    n1_filtered = Int[]
    n2_filtered = Int[]  
    n3_filtered = Int[]
    for x in n1
        x_coord_min = a1 * (x + minimum(basis[:, 1]))
        x_coord_max = a1 * (x + maximum(basis[:, 1]))
        if x_coord_max >= x_min && x_coord_min <= x_max
            push!(n1_filtered, x)
        end
    end
    for y in n2
        y_coord_min = a2 * (y + minimum(basis[:, 2]))
        y_coord_max = a2 * (y + maximum(basis[:, 2]))
        if y_coord_max >= y_min && y_coord_min <= y_max
            push!(n2_filtered, y)
        end
    end
    for z in n3
        z_coord_min = a3 * (z + minimum(basis[:, 3]))
        z_coord_max = a3 * (z + maximum(basis[:, 3]))
        if z_coord_max >= z_min && z_coord_min <= z_max
            push!(n3_filtered, z)
        end
    end
    # Generate atoms with optimized coordinate calculation
    for x in n1_filtered, y in n2_filtered, z in n3_filtered, i in 1:length(basisTypes)
        # Direct scalar multiplication instead of matrix operation
        coordinate = [a1 * (x + basis[i, 1]),
                     a2 * (y + basis[i, 2]), 
                     a3 * (z + basis[i, 3])]
        # Check if atom is within cell boundaries
        if (coordinate[1] >= x_min && coordinate[1] <= x_max &&
            coordinate[2] >= y_min && coordinate[2] <= y_max &&
            coordinate[3] >= z_min && coordinate[3] <= z_max)
            atom = Atom(basisTypes[i], coordinate, parameters)
            atom.latticeCoordinate .= atom.coordinate
            atom.cellIndex = cell.index
            atom.index = 0  # temporary value
            atom.isNewlyLoaded = true
            push!(atoms, atom)
        end
    end
    return atoms 
end

function ComputeLatticeAtoms_General!(cell::Cell, simulator::Simulator)
    parameters = simulator.parameters
    primaryVectors = parameters.primaryVectors
    latticeRanges = parameters.latticeRanges
    basisTypes = parameters.basisTypes
    basis = parameters.basis
    if !cell.isSavedAtomRange
        # Original calculation for non-orthogonal case
        primaryVectors_INV = parameters.primaryVectors_INV
        vertexMatrix = cell.vertexMatrix
        nfrac = vertexMatrix * primaryVectors_INV
        
        nmin = [floor(Int, minimum(nfrac[:, 1]) - 1), 
                floor(Int, minimum(nfrac[:, 2]) - 1), 
                floor(Int, minimum(nfrac[:, 3]) - 1)]
        nmax = [ceil(Int, maximum(nfrac[:, 1]) + 1),
                ceil(Int, maximum(nfrac[:, 2]) + 1),
                ceil(Int, maximum(nfrac[:, 3]) + 1)]
        
        n1 = max(nmin[1],latticeRanges[1,1]):min(nmax[1],latticeRanges[1,2])
        n2 = max(nmin[2],latticeRanges[2,1]):min(nmax[2],latticeRanges[2,2])
        n3 = max(nmin[3],latticeRanges[3,1]):min(nmax[3],latticeRanges[3,2])
        cell.atomRange = [n1, n2, n3]
        cell.isSavedAtomRange = true
    else
        n1, n2, n3 = cell.atomRange
    end
    # Get cell ranges for coordinate filtering
    ranges = cell.ranges
    x_min, x_max = ranges[1, 1], ranges[1, 2]
    y_min, y_max = ranges[2, 1], ranges[2, 2] 
    z_min, z_max = ranges[3, 1], ranges[3, 2]
    atoms = Vector{Atom}()
    # Original generation for non-orthogonal case
    for x in n1, y in n2, z in n3, i in 1:length(basisTypes)
        coordinate = primaryVectors' * (Float64[x, y, z] + basis[i, :])
        
        # Check if atom is within cell boundaries
        if (coordinate[1] >= x_min && coordinate[1] <= x_max &&
            coordinate[2] >= y_min && coordinate[2] <= y_max &&
            coordinate[3] >= z_min && coordinate[3] <= z_max)
            atom = Atom(basisTypes[i], coordinate, parameters)
            atom.latticeCoordinate .= atom.coordinate
            atom.cellIndex = cell.index
            atom.index = 0  # temporary value
            atom.isNewlyLoaded = true
            push!(atoms, atom)
        end
    end
    return atoms
end

function LoadCellAtoms!(cell::Cell, simulator::Simulator)
    if !cell.isLoaded
        if simulator.parameters.isPrimaryVectorOrthogonal
            atoms = ComputeLatticeAtoms_Orthogonal!(cell, simulator)
        else
            atoms = ComputeLatticeAtoms_General!(cell, simulator)
        end
        for vacancy in cell.vacancies
            for na in length(atoms):-1:1
                if ComputeDistance_squared(atoms[na].coordinate, vacancy.coordinate, (Int8(0), Int8(0), Int8(0)), simulator.box) < 1E-10
                    deleteat!(atoms, na)
                    break
                end
            end
        end
        # Apply perturbations
        for atom in atoms
            Pertubation!(atom, simulator)
        end
        cell.atomicDensity = length(atoms) / simulator.grid.cellVolume
        cell.latticeAtoms = atoms
        cell.isLoaded = true
    end
end

function GetTargetsFromNeighbor_dynamicLoad(atom::Atom, cell::Cell, filterIndexes::Vector{Int64}, simulator::Simulator)
    grid = simulator.grid
    box = simulator.box
    targets = Vector{Atom}()
    pMax = simulator.parameters.pMax
    nthreads = Threads.nthreads()
    cands_tls = [Atom[] for _ in 1:nthreads]
    neighborCellsInfo = cell.neighborCellsInfo
    AlreadyLoadedFlags = [true for _ in 1:27]
    infiniteFlag_tls = [true for _ in 1:27]

    if Threads.nthreads() > 1 
        for n in 1:length(neighborCellsInfo)
            neighborCellInfo = neighborCellsInfo[n]
            index = neighborCellInfo.index
            GetCell(grid, index)  # preload the cell to avoid race condition
        end
    end
    @threads :static for n in 1:length(neighborCellsInfo)
        neighborCellInfo = neighborCellsInfo[n]
        buf = cands_tls[Threads.threadid()]
        index = neighborCellInfo.index
        neighborCell = GetCell(grid, index)
        if neighborCell.isExplored
            continue
        end
        if !neighborCell.isLoaded
            AlreadyLoadedFlags[n] = false
        end
        LoadCellAtoms!(neighborCell, simulator)
        neighborCell.isExplored = true
        infiniteFlag_tls[n] = false
        for neighborAtom in [neighborCell.atoms; neighborCell.latticeAtoms]
            if neighborAtom.index == atom.index || neighborAtom.index in filterIndexes    
                continue
            end
            if ComputeVDistance(atom, neighborAtom, neighborCellInfo.cross, box) > 0 
                p = ComputeP!(atom, neighborAtom, neighborCellInfo.cross, box, pMax)
                if p >= pMax
                    continue
                end
                push!(buf, neighborAtom)
            end
        end
    end

    candidateTargets = reduce(vcat, cands_tls)
    infiniteFlag = reduce(&, infiniteFlag_tls)


    for (neighborCellInfo, flag) in zip(neighborCellsInfo, AlreadyLoadedFlags)
        idx = neighborCellInfo.index
        cell = GetCell(simulator.grid, idx)
        push!(simulator.exploredCells, cell)
        if flag == false
            push!(simulator.loadedCells, cell)
            for atom in cell.latticeAtoms
                simulator.minLatticeAtomID -= 1
                atom.index = simulator.minLatticeAtomID
            end
        end
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
        matchFlag = true
        for target in targets   
            if !SimultaneousCriteria(atom, candidateTarget, target, simulator)
                matchFlag = false
                break
            end
        end
        if matchFlag
            push!(targets, candidateTarget)
        end
    end    
    return (targets, infiniteFlag)
end

function Collision_dynamicLoad!(atom_p::Atom, atoms_t::Vector{Atom}, simulator::Simulator)
    N_t = length(atoms_t)
    grid = simulator.grid
    tanφList = Vector{Float64}(undef, N_t)
    tanψList = Vector{Float64}(undef, N_t)
    E_tList = Vector{Float64}(undef, N_t)
    x_pList = Vector{Float64}(undef, N_t)
    x_tList = Vector{Float64}(undef, N_t)
    Q_locList = Vector{Float64}(undef, N_t)
    pL = 0.0
    for atom_t in atoms_t
        l = atom_t.pL
        pL += l
    end
    pL /= N_t   
    if atom_p.numberOfEmptyCells > 1
        # This is an approximation, the pL is not accurate, but for most situations in bulk simulation, numberOfEmptyCells is 0.
        pL *= 1 / atom_p.numberOfEmptyCells
    end
    atom_t = atoms_t[1]
    N = GetCell(grid, atom_t.cellIndex).atomicDensity
    Q_nl_v = Q_nl(atom_p.energy, atom_p.mass, atom_t.mass, atom_p.type, atom_t.type,
                         pL, N, simulator.constantsByType)
    atom_p.energy -= Q_nl_v
    if atom_p.energy < 0.1 && atom_p.energy + Q_nl_v >= 0.1
        atom_p.energy = 0.11
    end
    for (i, atom_t) in enumerate(atoms_t)
        p = atom_t.pValue
        N = GetCell(grid, atom_t.cellIndex).atomicDensity 
        tanφList[i], tanψList[i], E_tList[i], x_pList[i], x_tList[i], Q_locList[i] = CollisionParams(
            atom_p.energy, atom_p.mass, atom_t.mass, atom_p.type, atom_t.type, p, simulator.constantsByType,
            simulator.θFunctions[[atom_p.type, atom_t.type]], simulator.τFunctions[[atom_p.type, atom_t.type]])
    end
    sumE_t = sum(E_tList)
    sumQ_loc = sum(Q_locList)
    η = N_t * atom_p.energy / (N_t * atom_p.energy + (N_t - 1) * (sumE_t + sumQ_loc))
    E_tList *= η   
    avePPoint = Vector{Float64}([0.0,0.0,0.0])
    momentum = Vector{Float64}([0.0,0.0,0.0])
    for (i, atom_t) in enumerate(atoms_t)
        if atom_t.pValue != 0
            velocityDirectionTmp = -atom_t.pVector / atom_t.pValue * tanψList[i] + atom_p.velocityDirection
        else
            velocityDirectionTmp = atom_p.velocityDirection
        end   
        SetVelocityDirection!(atom_t, velocityDirectionTmp)
        if E_tList[i] > GetDTE(atom_t, simulator) && E_tList[i] - GetBDE(atom_t, simulator) > 0.1
            LeaveLatticePoint_dynamicLoad!(atom_t, simulator)
            SetEnergy!(atom_t, E_tList[i] - GetBDE(atom_t, simulator))
            tCoordinate = atom_t.coordinate + x_tList[i] * η * atom_p.velocityDirection
            DisplaceAtom!(atom_t, tCoordinate, simulator)  
        else
            SetEnergy!(atom_t, 0.0)
        end
        avePPoint += atom_t.pPoint
        momentum += sqrt(2 * atom_t.mass * E_tList[i]) * atom_t.velocityDirection
    end

    # Update atom_p
    avePPoint /= N_t
    x_p = η * sum(x_pList) / N_t  # important modification
    pCoordinate = avePPoint - x_p * atom_p.velocityDirection
    #DisplaceAtom!(atom_p, avePPoint, simulator)
    DisplaceAtom!(atom_p, pCoordinate, simulator)
    velocity = (sqrt(2 * atom_p.mass * atom_p.energy) * atom_p.velocityDirection - momentum)  / atom_p.mass
    SetVelocityDirection!(atom_p, velocity)
    SetEnergy!(atom_p, atom_p.energy - (sumE_t + sumQ_loc) * η)
end 

function DumpInCascade_dynamicLoad(simulator::Simulator)
    if simulator.parameters.isDumpInCascade
        @dump "Cascade_$(simulator.nCascade).dump" [simulator.atoms; simulator.vacancies] []
        #Dump_dynamicLoad(simulator, "$(simulator.parameters.dumpFolder)/Cascade_$(simulator.nCascade).dump", simulator.nCollisionEvent, type)
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
                        DisplaceAtom!(target, target.coordinate, simulator)
                        push!(nextPAtoms, target)
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
        for cell in simulator.loadedCells
            cell.isLoaded = false
            empty!(cell.latticeAtoms) 
        end
        empty!(simulator.loadedCells)
        simulator.minLatticeAtomID = 0
    end
end

function delete_dynamicLoad!(simulator::Simulator, atom::Atom; isDeleteVacancy::Bool = false)
    cell = GetCell(simulator.grid, atom.cellIndex)
    if !isDeleteVacancy
        deleteat!(cell.atoms, findfirst(a -> a.index == atom.index, cell.atoms))
        simulator.numberOfAtoms -= 1
    else
        deleteat!(cell.vacancies, findfirst(v -> v.index == atom.index, cell.vacancies))
        if cell.isLoaded
            latom = CopyAtom(atom, simulator)
            simulator.minLatticeAtomID -= 1
            latom.index = simulator.minLatticeAtomID
            push!(cell.latticeAtoms, latom)
        end 
        simulator.numberOfVacancies -= 1
    end
    atom.isAlive = false 
end

function Stop_dynamicLoad!(atom::Atom, simulator::Simulator)
    grid = simulator.grid
    cell = GetCell(grid, atom.cellIndex)
    nearestVacancyDistance_squared = Inf
    isExist = false
    nearestVacancy = Atom(1, [0.0,0.0,0.0], simulator.parameters)  # declare variable to store the nearest vacancy
    if ! cell.isPushedNeighbor
        SetCellNeighborInfo!(cell, grid)
        cell.isPushedNeighbor = true
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
    if isExist 
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
    vacancy = Atom(atom.type, atom.latticeCoordinate[:], simulator.parameters)
    vacancy.cellIndex = atom.cellIndex
    vacancy.type += length(keys(simulator.parameters.typeDict))
    return vacancy
end



function ShotTarget_dynamicLoad(atom::Atom, filterIndexes::Vector{Int64}, simulator::Simulator)
    grid = simulator.grid
    periodic = simulator.parameters.periodic    
    cell = GetCell(grid, atom.cellIndex)
    atom.numberOfEmptyCells = 0
    while true
        if ! cell.isPushedNeighbor
            SetCellNeighborInfo!(cell, grid)
            cell.isPushedNeighbor = true
        end
        targets, isInfinity = GetTargetsFromNeighbor_dynamicLoad(atom, cell, filterIndexes, simulator)
        if length(targets) > 0
            for cell in simulator.exploredCells
                cell.isExplored = false
            end
            empty!(simulator.exploredCells)
            return targets, true
        else
            atom.numberOfEmptyCells += 1
            dimension, direction = AtomOutFaceDimension(atom, cell)
            neighborIndex = Vector{Int8}([0,0,0])
            neighborIndex[dimension] = direction == 1 ? Int8(-1) : Int8(1)
            neighborIndex .+= 2
            neighborInfo = cell.neighborCellsInfo[neighborIndex...]
            crossFlag = neighborInfo.cross
            if crossFlag[dimension] != 0 && periodic[dimension]
                atom.coordinate[dimension] -= crossFlag[dimension] * simulator.box.vectors[dimension, dimension]
            end
            if (neighborInfo.cross[dimension] != 0 && !periodic[dimension]) || isInfinity
                for cell in simulator.exploredCells
                    cell.isExplored = false
                end
                atom.numberOfEmptyCells = 0
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
            if cell.isLoaded && atom.type > length(keys(simulator.parameters.typeDict)) 
                latticeAtom = Atom(atom.type-length(keys(simulator.parameters.typeDict)), atom.coordinate, parameters)
                latticeAtom.latticeCoordinate .= atom.coordinate
                Pertubation!(latticeAtom, simulator)
                latticeAtom.cellIndex = cell.index
                simulator.minLatticeAtomID -= 1
                latticeAtom.index = simulator.minLatticeAtomID
                latticeAtom.isNewlyLoaded = true
                push!(cell.latticeAtoms, latticeAtom)
            end
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
