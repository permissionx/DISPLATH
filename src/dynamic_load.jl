function LoadCell(cell::GridCell, simulator::Simulator)
    if !cell.isLoaded
        #println(cell.index)
        # only for single type of target material
        parameters = simulator.parameters
        primaryVectors = parameters.primaryVectors
        latticeRanges = parameters.latticeRanges
        basisTypes = parameters.basisTypes
        basis = parameters.basis

        if !cell.isSavedAtomRange
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
        
        atoms = Vector{Atom}()  # Use dynamic array instead of pre-allocation
        
        # Generate atoms and filter by cell boundaries
        for x in n1 
            for y in n2
                for z in n3
                    for i in  basisTypes
                        coordinate = primaryVectors' * (Float64[x, y, z] + basis[i, :])
                        
                        # Check if atom is within cell boundaries
                        if (coordinate[1] >= x_min && coordinate[1] <= x_max &&
                            coordinate[2] >= y_min && coordinate[2] <= y_max &&
                            coordinate[3] >= z_min && coordinate[3] <= z_max)
                            
                            atom = Atom(basisTypes[i], coordinate, parameters)
                            atom.cellIndex = cell.index
                            simulator.minLatticeAtomID -= 1
                            atom.index = simulator.minLatticeAtomID
                            push!(atoms, atom)
                        end
                    end
                end
            end
        end
        
        for vacancy in cell.vacancies
            for na in length(atoms):-1:1
                if ComputeDistance_squared(atoms[na].coordinate, vacancy.coordinate, Int8[0,0,0] ,simulator.box) < 1E-10
                    deleteat!(atoms, na)
                    break
                end
            end
        end
        cell.latticeAtoms = atoms
        cell.isLoaded = true
        push!(simulator.loadedCells, cell)
    end
    cell.allAtoms = [cell.latticeAtoms; cell.atoms]
    #println(length(cell.allAtoms))
end

function GetTargetsFromNeighbor_dynamicLoad(atom::Atom, gridCell::GridCell, simulator::Simulator)
    # todo: This is not prices for atom may be added to targets with lower height
    cellGrid = simulator.cellGrid
    box = simulator.box
    targets = Vector{Atom}()
    infiniteFlag = true
    for (_, neighborCellInfo) in gridCell.neighborCellsInfo
        index = neighborCellInfo.index
        neighborCell = cellGrid.cells[index[1], index[2], index[3]]
        if neighborCell.isExplored
            continue
        end
        LoadCell(neighborCell, simulator)
        infiniteFlag = false
        for neighborAtom in neighborCell.allAtoms
            if ComputeVDistance(atom, neighborAtom, neighborCellInfo.cross, box) > 0 && neighborAtom.index != atom.index
                originalCoordinate = copy(neighborAtom.coordinate)
                Pertubation!(neighborAtom, simulator)
                ComputeP!(atom, neighborAtom, neighborCellInfo.cross, box)
                if neighborAtom.pValue[atom.index] >= simulator.parameters.pMax
                    DeleteP!(neighborAtom, atom.index)
                    continue
                end
                matchFlag = true
                for target in targets
                    if !SimultaneousCriteria(atom, neighborAtom, target, neighborCellInfo.cross, simulator)
                        matchFlag = false
                        DeleteP!(neighborAtom, atom.index)
                        neighborAtom.coordinate = originalCoordinate[:]
                        break
                    end
                end
                if matchFlag
                    push!(targets, neighborAtom)
                end
            end
        end
        neighborCell.isExplored = true
        push!(simulator.exploredCells, neighborCell)
    end
    if infiniteFlag
        log("Infinitely fly atom in the $(simulator.nIrradiation)th irradiation:\n$(atom)\n")
    end
    return (targets, infiniteFlag)
end


function Collision_dynamicLoad!(atom_p::Atom, atoms_t::Vector{Atom}, simulator::Simulator, nStep::Int64)
    N_t = length(atoms_t)
    cellGrid = simulator.cellGrid
    tanφList = Vector{Float64}(undef, N_t)
    tanψList = Vector{Float64}(undef, N_t)
    E_tList = Vector{Float64}(undef, N_t)
    x_pList = Vector{Float64}(undef, N_t)
    x_tList = Vector{Float64}(undef, N_t)
    QList = Vector{Float64}(undef, N_t)
    pL = 0.0
    for atom_t in atoms_t
        l = atom_t.pL[atom_p.index]
        #if l > simulator.parameters.pLMax 
        #    l = simulator.parameters.pLMax
        #end
        pL += l
    end
    pL /= N_t   
    for (i, atom_t) in enumerate(atoms_t)
        p = atom_t.pValue[atom_p.index]
        N = cellGrid.cells[atom_t.cellIndex[1], atom_t.cellIndex[2], atom_t.cellIndex[3]].atomicDensity 
        tanφList[i], tanψList[i], E_tList[i], x_pList[i], x_tList[i], QList[i] = CollisionParams(
            atom_p.energy, atom_p.mass, atom_t.mass, atom_p.type, atom_t.type, p, pL, N, simulator.constantsByType,
            simulator.θFunctions[[atom_p.type, atom_t.type]], simulator.τFunctions[[atom_p.type, atom_t.type]])
    end
    sumE_t = sum(E_tList)
    η = N_t * atom_p.energy / (N_t * atom_p.energy + (N_t - 1) * sumE_t)
    E_tList *= η
    # Update atoms_t (target atoms)         
    avePPoint = Vector{Float64}([0.0,0.0,0.0])
    momentum = Vector{Float64}([0.0,0.0,0.0])
    for (i, atom_t) in enumerate(atoms_t)
        if atom_t.pValue[atom_p.index] != 0
            velocityDirectionTmp = -atom_t.pVector[atom_p.index] / atom_t.pValue[atom_p.index] * tanψList[i] + atom_p.velocityDirection
        else
            velocityDirectionTmp = atom_p.velocityDirection
        end   
        SetVelocityDirection!(atom_t, velocityDirectionTmp)
        if E_tList[i] > GetDTE(atom_t, simulator) && E_tList[i] - GetBDE(atom_t, simulator) > 0.0
            LeaveLatticePoint_dynamicLoad!(atom_t, simulator)
            SetEnergy!(atom_t, E_tList[i] - GetBDE(atom_t, simulator))
            tCoordinate = atom_t.coordinate + x_tList[i] * η * atom_p.velocityDirection
            DisplaceAtom_dynamicLoad!(atom_t, tCoordinate, simulator)  
            SetEnergy!(atom_t, E_tList[i])
        end
        avePPoint += atom_t.pPoint[atom_p.index]
        momentum += sqrt(2 * atom_t.mass * E_tList[i]) * atom_t.velocityDirection
    end

    # Update atom_p
    avePPoint /= N_t
    x_p = η * sum(x_pList)
    pCoordinate = avePPoint - x_p * atom_p.velocityDirection
    DisplaceAtom_dynamicLoad!(atom_p, pCoordinate, simulator)
    velocity = (sqrt(2 * atom_p.mass * atom_p.energy) * atom_p.velocityDirection - momentum)  / atom_p.mass
    SetVelocityDirection!(atom_p, velocity)
    
    SetEnergy!(atom_p, atom_p.energy - (sumE_t + sum(QList)) * η)
end 

function DisplaceAtom_dynamicLoad!(atom::Atom, coordinate::Vector{Float64}, simulator::Simulator)
    cell = simulator.cellGrid.cells[atom.cellIndex[1], atom.cellIndex[2], atom.cellIndex[3]]
    d = cell.atomicDensity
    DisplaceAtom!(atom, coordinate, simulator)
    cell.atomicDensity = d
end

function Cascade_dynamicLoad!(atom_p::Atom, simulator::Simulator)
    pAtoms = Vector{Atom}([atom_p])

    while true
        targetsList = Vector{Vector{Atom}}()
        for atom in pAtoms
            targets = ShotTarget_dynamicLoad(atom, simulator)
            for pAtom in pAtoms
                filter!(t->t.index != pAtom.index, targets)
            end
            push!(targetsList, targets)
        end
        UniqueTargets!(targetsList, pAtoms)
        nextPAtoms = Vector{Atom}()
        for (pAtom, targets) in zip(pAtoms, targetsList)
            if length(targets) == 0.0
                pAtom.lastTargets = Vector{Int64}()
                delete_dynamicLoad!(simulator, pAtom)
                continue
            end
            pAtom.lastTargets = [t.index for t in targets]
            Collision_dynamicLoad!(pAtom, targets, simulator, simulator.nDisplaceStep)
            for target in targets
                if target.energy > 0.0   
                    push!(nextPAtoms, target)
                end
            end
            if pAtom.energy > simulator.parameters.stopEnergy   
                push!(nextPAtoms, pAtom)
            else
                pAtom.lastTargets = Vector{Int64}()
                Stop_dynamicLoad!(pAtom, simulator)
            end
        end
        for targets in targetsList
            for target in targets
                EmptyP!(target)
            end
        end 
        simulator.nDisplaceStep += 1
        if simulator.parameters.isDumpInCascade
            Dump_dynamicLoad(simulator, "Cascade.dump", simulator.nDisplaceStep, true)
        end
        if simulator.parameters.isLog
            println("Collision times: ", simulator.nDisplaceStep)
        end
        if length(nextPAtoms) > 0
            pAtoms = nextPAtoms
        else
            break
        end
    end
    for cell in simulator.loadedCells
        cell.isLoaded = false
        cell.allAtoms = Vector{Atom}()
    end
    empty!(simulator.loadedCells)
    simulator.minLatticeAtomID = 0
end

function delete_dynamicLoad!(simulator::Simulator, atom::Atom; isDeleteVacancy::Bool = false)
    cell = simulator.cellGrid.cells[atom.cellIndex[1], atom.cellIndex[2], atom.cellIndex[3]]
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
    cellGrid = simulator.cellGrid
    cell = cellGrid.cells[atom.cellIndex[1], atom.cellIndex[2], atom.cellIndex[3]]
    nearestVacancyDistance_squared = Inf
    isExist = false
    nearestVacancy = Atom(1, [0.0,0.0,0.0], simulator.parameters)  # declare variable to store the nearest vacancy
    if ! cell.isPushedNeighbor
        IterPushCellNeighbors!(cellGrid, cell)
        cell.isPushedNeighbor = true
    end
    for (_ , neighborCellInfo) in cell.neighborCellsInfo
        index = neighborCellInfo.index
        cross = neighborCellInfo.cross
        neighborCell = simulator.cellGrid.cells[index[1], index[2], index[3]]
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
        if atom.type == nearestVacancy.type
            delete_dynamicLoad!(simulator, atom)
            delete_dynamicLoad!(simulator, nearestVacancy, isDeleteVacancy = true)
        else
            atom.coordinate = nearestVacancy.coordinate[:]
        end
    end
end




function LeaveLatticePoint_dynamicLoad!(atom::Atom, simulator::Simulator; isUpdateEnv::Bool = true)
    if atom.isNewlyLoaded
        cell = simulator.cellGrid.cells[atom.cellIndex[1], atom.cellIndex[2], atom.cellIndex[3]]
        atom.isNewlyLoaded = false
        vacancy = copyAtom(atom, simulator)
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

function copyAtom(atom::Atom, simulator::Simulator)
    newAtom = Atom(atom.type, atom.coordinate, simulator.parameters)
    newAtom.cellIndex = atom.cellIndex[:]
    return newAtom
end



function ShotTarget_dynamicLoad(atom::Atom, simulator::Simulator)
    cellGrid = simulator.cellGrid
    periodic = simulator.parameters.periodic    
    cell = cellGrid.cells[atom.cellIndex[1], atom.cellIndex[2], atom.cellIndex[3]]
    #accCrossFlag = Vector{Int64}([0,0,0])
    while true
        if ! cell.isPushedNeighbor
            IterPushCellNeighbors!(cellGrid, cell)
            cell.isPushedNeighbor = true
        end
        (targets, isInfinity) = GetTargetsFromNeighbor_dynamicLoad(atom, cell, simulator)
        # delete repeated targets in lastTargets
        filter!(t->!(t.index in atom.lastTargets), targets)
        if length(targets) > 0
            for cell in simulator.exploredCells
                cell.isExplored = false
            end
            empty!(simulator.exploredCells)
            return targets
        else
            dimension, direction = AtomOutFaceDimension(atom, cell)
            neighborIndex = Vector{Int8}([0,0,0])
            neighborIndex[dimension] = direction == 1 ? Int8(-1) : Int8(1)
            neighborInfo = cell.neighborCellsInfo[neighborIndex]
            crossFlag = neighborInfo.cross
            for d in 1:3
                if crossFlag[d] != 0
                    atom.coordinate[d] -= crossFlag[d] * simulator.box.vectors[d, d]
                end
            end 
            if (neighborInfo.cross[dimension] != 0 && !periodic[dimension]) || isInfinity
                for cell in simulator.exploredCells
                    cell.isExplored = false
                end
                empty!(simulator.exploredCells)
                return Vector{Atom}() # means find nothing  
            end 
            index = neighborInfo.index
            cell = cellGrid.cells[index[1], index[2], index[3]]
        end
    end
end


function Dump_dynamicLoad(simulator::Simulator, fileName::String, step::Int64, isAppend::Bool=false, isDebug::Bool=false)
    if !simulator.parameters.isOrthogonal        
        error("The box is not orthogonal, please use the orthogonal box.")
    end
    write_flag = isAppend ? "a" : "w"
    open(fileName, write_flag) do file
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
                write(file, "$(atom.index+100000) $(atom.type+5) \
                $(atom.coordinate[1]) $(atom.coordinate[2]) $(atom.coordinate[3]) 0.0\n")
            end
        end
    end
end