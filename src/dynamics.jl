using .BCA
using Printf

function ShotTarget(atom::Atom, filterIndexes::Vector{Int64}, simulator::Simulator)
    grid = simulator.grid
    periodic = simulator.parameters.periodic    
    cell = GetCell(grid, atom.cellIndex)
    atom.emptyPath = 0.0
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
            dimension, direction,t = AtomOutFaceDimension(atom, cell)
            atom.emptyPath += t
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
                atom.emptyPath = 0.0
                empty!(simulator.exploredCells)
                return Vector{Atom}(), false # means find nothing  
            end 
            index = neighborInfo.index
            cell = GetCell(grid, index)
        end
    end
end




function AtomOutFaceDimension(atom::Atom, cell::Cell)
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
            return d, rangeIndex, t
        end
    end
    error("Out face not found\n ########Atom#######\n $(atom) \n 
                                ########cell#######\n $(cell.ranges) \n $(cell.index)\n")
end


function GetTargetsFromNeighbor(atom::Atom, cell::Cell, filterIndexes::Vector{Int64}, simulator::Simulator)
    grid = simulator.grid
    box = simulator.box
    targets = Vector{Atom}()
    infiniteFlag = true
    candidateTargets = Vector{Atom}()
    pMax = simulator.parameters.pMax
    for neighborCellInfo in cell.neighborCellsInfo
        index = neighborCellInfo.index
        neighborCell = GetCell(grid, index)
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
            if ComputeVDistance(atom, neighborAtom, neighborCellInfo.cross, box) > 0 
                p = ComputeP!(atom, neighborAtom, neighborCellInfo.cross, box)
                if p >= pMax
                    continue
                end
                push!(candidateTargets, neighborAtom)
            end
        end
    end
    if infiniteFlag
        @record "log" "Infinitely fly atom in the $(simulator.nCascade)th irradiation:\n$(atom)"
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
            if !SimultaneousCriteria(candidateTarget, target, simulator)
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


function Collision!(atom_p::Atom, atoms_t::Vector{Atom}, simulator::Simulator)
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
    pL -= atom_p.emptyPath
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
            SetEnergy!(atom_t, E_tList[i] - GetBDE(atom_t, simulator))
            tCoordinate = atom_t.coordinate + x_tList[i] * η * atom_p.velocityDirection
            DisplaceAtom!(atom_t, tCoordinate, simulator)  
            if atom_t.latticePointIndex != -1
                LeaveLatticePoint!(atom_t, simulator)
            end     
        else 
            SetEnergy!(atom_t, 0.0)
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
    if !IS_DYNAMIC_LOAD
        Cascade_staticLoad!(atom_p, simulator)
    else
        Cascade_dynamicLoad!(atom_p, simulator)
    end
end


function Cascade_staticLoad!(atom_p::Atom, simulator::Simulator)
    pAtoms = Vector{Atom}([atom_p])
    pAtomsIndex = [a.index for a in pAtoms]
    parameters = simulator.parameters
    simulator.nCollisionEvent = 0
    simulator.nCascade += 1
    DumpInCascade(simulator)
    while true
        simulator.nCollisionEvent += 1
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
            if length(targets) > 0
                pAtom.lastTargets = [t.index for t in targets]
                Collision!(pAtom, targets, simulator)
                for target in targets
                    if target.energy > 0.0   
                        #DisplaceAtom!(target, target.coordinate, simulator)
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
        DumpInCascade(simulator)
        if length(nextPAtoms) > 0
            pAtoms = nextPAtoms
            sort!(pAtoms, by = a -> a.energy, rev = true)
            pAtomsIndex = [a.index for a in pAtoms]
        else
            break
        end
    end
end

function DumpInCascade(simulator::Simulator)
    if simulator.parameters.isDumpInCascade
        @dump "Cascade_$(simulator.nCascade).dump" simulator.atoms ["vx", "vy", "vz", "e"]
    end
end

