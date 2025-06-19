using .BCA
using Printf

function ShotTarget(atom::Atom, filterIndexes::Vector{Int64}, simulator::Simulator)
    cellGrid = simulator.cellGrid
    periodic = simulator.parameters.periodic    
    cell = cellGrid.cells[atom.cellIndex...]
    while true
        (targets, isInfinity) = GetTargetsFromNeighbor(atom, cell, filterIndexes, simulator)
        # delete repeated targets in lastTargets
        if length(targets) > 0
            for cell in simulator.exploredCells
                cell.isExplored = false
            end
            empty!(simulator.exploredCells)
            return targets, true
        else
            dimension, direction = AtomOutFaceDimension(atom, cell)
            neighborIndex = Vector{Int8}([0,0,0])
            neighborIndex[dimension] = direction == 1 ? Int8(-1) : Int8(1)
            neighborIndex .+= 2
            neighborInfo = cell.neighborCellsInfo[neighborIndex...]
            crossFlag = neighborInfo.cross
            if crossFlag[dimension] != 0 && periodic[dimension]
                atom.coordinate[dimension] -= crossFlag[dimension] * simulator.box.vectors[dimension, dimension]
            end
            if (crossFlag[dimension] != 0 && !periodic[dimension]) || isInfinity
                for cell in simulator.exploredCells
                    cell.isExplored = false
                end
                empty!(simulator.exploredCells)
                return Vector{Atom}(), false # means find nothing  
            end 
            index = neighborInfo.index
            cell = cellGrid.cells[index...]
        end
    end
end


function AtomOutFaceDimension(atom::Atom, cell::GridCell)
    coordinate = atom.coordinate

    for d in 1:3
        if atom.velocityDirection[d] >= 0
            rangeIndex = 2
        else
            rangeIndex = 1
        end
        faceCoordinate = cell.ranges[d, rangeIndex]
        t = (faceCoordinate - coordinate[d]) / atom.velocityDirection[d]
        elseDs = [ed for ed in 1:3 if ed != d]
        allInRange = true
        for elseD in elseDs
            crossCoord = coordinate[elseD] + atom.velocityDirection[elseD] * t
            if !(cell.ranges[elseD, 1] <= crossCoord <= cell.ranges[elseD, 2])
                allInRange = false
                break
            end
        end
        if allInRange
            return d, rangeIndex
        end
    end
    error("Out face not found\n ########Atom#######\n $(atom) \n 
                                ########cell#######\n $(cell.ranges) \n $(cell.index)\n")
end


function Collision!(atom_p::Atom, atoms_t::Vector{Atom}, simulator::Simulator)
    N_t = length(atoms_t)
    cellGrid = simulator.cellGrid
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
        pL *= (atom_p.numberOfEmptyCells - 1) / atom_p.numberOfEmptyCells
    end
    atom_t = atoms_t[1]
    N = cellGrid.cells[atom_t.cellIndex[1], atom_t.cellIndex[2], atom_t.cellIndex[3]].atomicDensity
    Q_nl_v = Q_nl(atom_p.energy, atom_p.mass, atom_t.mass, atom_p.type, atom_t.type,
                         pL, N, simulator.constantsByType)  
    atom_p.energy -= Q_nl_v
    for (i, atom_t) in enumerate(atoms_t)
        p = atom_t.pValue
        N = cellGrid.cells[atom_t.cellIndex[1], atom_t.cellIndex[2], atom_t.cellIndex[3]].atomicDensity 
        tanφList[i], tanψList[i], E_tList[i], x_pList[i], x_tList[i], Q_locList[i] = CollisionParams(
            atom_p.energy, atom_p.mass, atom_t.mass, atom_p.type, atom_t.type, p, simulator.constantsByType,
            simulator.θFunctions[[atom_p.type, atom_t.type]], simulator.τFunctions[[atom_p.type, atom_t.type]])
        if atom_p.type == 3
            Log("$(E_tList[i])\n", simulator, fileName="Et.csv")
        end
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
            SetEnergy!(atom_t, E_tList[i] - GetBDE(atom_t, simulator))
            tCoordinate = atom_t.coordinate + x_tList[i] * η * atom_p.velocityDirection
            DisplaceAtom!(atom_t, tCoordinate, simulator)  
            #SetEnergy!(atom_t, E_tList[i])
            if atom_t.latticePointIndex != -1
                LeaveLatticePoint!(atom_t, simulator)
            end     
        else 
            SetEnergy!(atom_t, 0.0)
            if atom_t.latticePointIndex != -1
                SetCoordinate!(atom_t, simulator.latticePoints[atom_t.latticePointIndex].coordinate)
            end
        end
        avePPoint += atom_t.pPoint
        momentum += sqrt(2 * atom_t.mass * E_tList[i]) * atom_t.velocityDirection
    end

    avePPoint /= N_t
    x_p = η * sum(x_pList) / N_t  # important modification
    pCoordinate = avePPoint - x_p * atom_p.velocityDirection
    DisplaceAtom!(atom_p, pCoordinate, simulator)
    velocity = (sqrt(2 * atom_p.mass * atom_p.energy) * atom_p.velocityDirection - momentum)  / atom_p.mass
    SetVelocityDirection!(atom_p, velocity)
    
    SetEnergy!(atom_p, atom_p.energy - (sumE_t + sumQ_loc) * η)
end 

function Cascade!(atom_p::Atom, simulator::Simulator)
    pAtoms = Vector{Atom}([atom_p])
    pAtomsIndex = [a.index for a in pAtoms]
    parameters = simulator.parameters
    simulator.nCollisionEvent = 0
    DumpInCascade(simulator, "w")
    while true
        targetsList = Vector{Vector{Atom}}()
        deleteIndexes = Int64[]
        othersTargetIndexes = Int64[]
        for (na, pAtom) in enumerate(pAtoms)
            targets, isAlive = ShotTarget(pAtom, [pAtomsIndex; pAtom.lastTargets; othersTargetIndexes], simulator)
            if !isAlive
                empty!(pAtom.lastTargets)
                delete!(simulator, pAtom)
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
            if pAtom.type == 3
                Log("$(length(targets))\n", simulator, fileName="targetNumber.csv")
                for target in targets
                    Log("$(target.pValue[1])\n", simulator, fileName="pValue.csv")
                end
            end
            if length(targets) > 0
                pAtom.lastTargets = [t.index for t in targets]
                Collision!(pAtom, targets, simulator)
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
                    Stop!(pAtom, simulator)
                end
            else
                push!(nextPAtoms, pAtom)
            end
        end
        for targets in targetsList
            for target in targets
                EmptyP!(target)
            end
        end 
        simulator.nCollisionEvent += 1
        DumpInCascade(simulator, "a")
        if length(nextPAtoms) > 0
            pAtoms = nextPAtoms
            sort!(pAtoms, by = a -> a.energy, rev = true)
            pAtomsIndex = [a.index for a in pAtoms]
        else
            break
        end
    end
    simulator.nCascade += 1
end





function GetTargetsFromNeighbor(atom::Atom, gridCell::GridCell, filterIndexes::Vector{Int64}, simulator::Simulator)
    cellGrid = simulator.cellGrid
    box = simulator.box
    targets = Vector{Atom}()
    infiniteFlag = true
    candidateTargets = Vector{Atom}()
    pMax = simulator.pMax
    for neighborCellInfo in gridCell.neighborCellsInfo
        index = neighborCellInfo.index
        neighborCell = cellGrid.cells[index...]
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
                p = ComputeP!(atom, neighborAtom, neighborCellInfo.cross, box, pMax)
                if p >= pMax
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
        Log("Infinitely fly atom in the $(simulator.nCascade)th irradiation:\n$(atom)\n", simulator)
    end
    if isempty(candidateTargets)
        return (targets, infiniteFlag)
    end
    # Find target with minimum pL value using Julia's built-in findmin
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
