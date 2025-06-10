using .BCA
using Printf

function ShotTarget(atom::Atom, filterIndexes::Vector{Int64}, simulator::Simulator)
    cellGrid = simulator.cellGrid
    periodic = simulator.parameters.periodic    
    cell = cellGrid.cells[atom.cellIndex[1], atom.cellIndex[2], atom.cellIndex[3]]
    while true
        (targets, isInfinity) = GetTargetsFromNeighbor(atom, cell, filterIndexes, simulator)
        # delete repeated targets in lastTargets
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
            if crossFlag[dimension] != 0 && periodic[dimension]
                atom.coordinate[dimension] -= crossFlag[dimension] * simulator.box.vectors[dimension, dimension]
            end
            if (crossFlag[dimension] != 0 && !periodic[dimension]) || isInfinity
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
        # debug: turn off the inelastic collision
        #println("DEBUG:\n tanφ: $(tanφList[i])\n tanψ: $(tanψList[i])\n E_t: $(E_tList[i])\n x_p: $(x_pList[i])\n x_t: $(x_tList[i])\n Q: $(QList[i])")
    end
    sumE_t = sum(E_tList)
    sumQ = sum(QList)
    η = N_t * atom_p.energy / (N_t * atom_p.energy + (N_t - 1) * (sumE_t+sumQ))
    E_tList *= η
    # Update atoms_t (target atoms)         
    avePPoint = Vector{Float64}([0.0,0.0,0.0])
    momentum = Vector{Float64}([0.0,0.0,0.0])
    # for cross section calculations
    #global tid
    #global buf
    #if (! (tid in [atom_t.index for atom_t in atoms_t])) && simulator.nCollisionEvent == 1
    #    write(buf, "0.0\n")
    #end
    for (i, atom_t) in enumerate(atoms_t)
        if atom_t.pValue[atom_p.index] != 0
            velocityDirectionTmp = -atom_t.pVector[atom_p.index] / atom_t.pValue[atom_p.index] * tanψList[i] + atom_p.velocityDirection
        else
            velocityDirectionTmp = atom_p.velocityDirection
        end   
        SetVelocityDirection!(atom_t, velocityDirectionTmp)
        # for cross section calculations
        #if atom_t.index == tid && simulator.nCollisionEvent == 1
        #    @printf(buf, "%f\n", E_tList[i])
        #end
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
            # this is for undoing the pertubation 
            if atom_t.latticePointIndex != -1
                SetCoordinate!(atom_t, simulator.latticePoints[atom_t.latticePointIndex].coordinate)
            end
        end
        avePPoint += atom_t.pPoint[atom_p.index]
        momentum += sqrt(2 * atom_t.mass * E_tList[i]) * atom_t.velocityDirection
    end

    # Update atom_p
    # initMomentum = sqrt(2 * atom_p.mass * atom_p.energy) * atom_p.velocityDirection   # Check momentum conservation 

    avePPoint /= N_t
    x_p = η * sum(x_pList)
    pCoordinate = avePPoint - x_p * atom_p.velocityDirection
    DisplaceAtom!(atom_p, pCoordinate, simulator)
    velocity = (sqrt(2 * atom_p.mass * atom_p.energy) * atom_p.velocityDirection - momentum)  / atom_p.mass
    SetVelocityDirection!(atom_p, velocity)
    
    SetEnergy!(atom_p, atom_p.energy - (sumE_t + sum(QList)) * η)
    #if atom_p.type == 2 
    #    Log("$(sumE_t),$(sum(QList))\n", simulator, fileName="loss/$(simulator.nCascade)")
    #end




    # Check momentum conservation 
    #afterMomentum = sqrt(2 * atom_p.mass * atom_p.energy) * atom_p.velocityDirection
    #for i in 1:length(atoms_t)
    #    afterMomentum += sqrt(2 * atoms_t[i].mass * E_tList[i]) * atoms_t[i].velocityDirection  # debug
    #end
    #open("p.debug.log", "a") do f
    #    write(f, "$(simulator.nCascade),$(N_t),$(initMomentum[1]),$(initMomentum[2]),$(initMomentum[3]),$(afterMomentum[1]),$(afterMomentum[2]),$(afterMomentum[3])\n")
    #end
    #exit()
end 




function Cascade!(atom_p::Atom, simulator::Simulator)
    pAtoms = Vector{Atom}([atom_p])
    pAtomsIndex = [a.index for a in pAtoms]
    if simulator.parameters.isDumpInCascade
        DumpDefects(simulator, "Cascade_$(simulator.nCascade).dump", simulator.nCollisionEvent, false)
    end
    while true
        targetsList = Vector{Vector{Atom}}()
        deleteIndexes = Vector{Int64}()
        for (na, pAtom) in enumerate(pAtoms)
            targets = ShotTarget(pAtom, [pAtomsIndex; pAtom.lastTargets], simulator)
            if length(targets) == 0.0
                pAtom.lastTargets = Vector{Int64}()
                delete!(simulator, pAtom)
                push!(deleteIndexes, na)
                continue
            end
            push!(targetsList, targets)
        end
        deleteat!(pAtoms, deleteIndexes)
        pAtomsIndex = [a.index for a in pAtoms]
        UniqueTargets!(targetsList, pAtoms)
        nextPAtoms = Vector{Atom}()
        for (pAtom, targets) in zip(pAtoms, targetsList)
            if length(targets) > 0
                pAtom.lastTargets = [t.index for t in targets]
                Collision!(pAtom, targets, simulator)
                for target in targets
                    if target.energy > 0.0   
                        DisplaceAtom!(target, target.coordinate, simulator)
                        push!(nextPAtoms, target)
                    end
                end
                if pAtom.energy > simulator.parameters.stopEnergy   
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
        if simulator.parameters.isDumpInCascade
            DumpDefects(simulator, "Cascade_$(simulator.nCascade).dump", simulator.nCollisionEvent, true)
        end
        if length(nextPAtoms) > 0
            pAtoms = nextPAtoms
            pAtomsIndex = [a.index for a in pAtoms]
        else
            break
        end
    end
    simulator.nCascade += 1
end


function UniqueTargets!(targetsList::Vector{Vector{Atom}}, pAtoms::Vector{Atom})
    targetToListDict = Dict{Int64, Vector{Int64}}()
    for (i, targets) in enumerate(targetsList)
        for target in targets
            try 
                push!(targetToListDict[target.index], i)
            catch 
                targetToListDict[target.index] = Vector{Int64}([i])
            end
        end
    end
    for (targetIndex, targetsListIndex) in targetToListDict
        if length(targetsListIndex) > 1
            maxEnergy = -1.0
            maxArg = 0
            for index in targetsListIndex
                energy = pAtoms[index].energy
                if energy > maxEnergy
                    maxEnergy = energy
                    maxArg = index
                end
            end
            for index in targetsListIndex
                if index != maxArg 
                    filter!(t -> t.index != targetIndex, targetsList[index])
                end
            end
        end
    end
end

