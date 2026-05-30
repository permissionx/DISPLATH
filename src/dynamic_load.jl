using StaticArrays

function ShotTarget_dynamicLoad(atom::Atom, filterIndexes::Vector{Int64}, simulator::Simulator)
    grid = simulator.grid
    periodic = simulator.parameters.periodic    
    cell = GetCell(grid, atom.cellIndex, simulator)
    emptyPath = 0.0
    while true
        targets = GetTargetsFromNeighbor_dynamicLoad(atom, cell, filterIndexes, simulator)
        if length(targets) > 0
            return targets, true, emptyPath
        else
            dimension, direction, t = AtomOutFaceDimension(atom, cell)
            emptyPath = t
            neighborIndex = MVector{3,Int8}(2, 2, 2)  
            neighborIndex[dimension] = direction == 1 ? Int8(1) : Int8(3)
            neighborCellsInfo = GetNeighborCellsInfo!(cell, grid, simulator)
            neighborInfo = neighborCellsInfo[neighborIndex[1], neighborIndex[2], neighborIndex[3]]
            crossFlag = neighborInfo.cross
            for i in 1:3
                for j in 1:3
                    oppoDirection = direction == 1 ? Int8(3) : Int8(1)
                    if dimension == 1
                        deprecatedCellKey = (Int(oppoDirection), i, j)
                    elseif dimension == 2
                        deprecatedCellKey = (i, Int(oppoDirection), j)
                    else
                        deprecatedCellKey = (i, j, Int(oppoDirection))
                    end
                    cell = GetCell(grid, neighborCellsInfo[deprecatedCellKey...].index, simulator)
                    DeprecateCell!(cell, simulator)
                end
            end
            if crossFlag[dimension] != 0 && periodic[dimension]
                atom.coordinate[dimension] -= crossFlag[dimension] * simulator.box.vectors[dimension, dimension]
            end
            if (neighborInfo.cross[dimension] != 0 && !periodic[dimension]) || t >= simulator.parameters.infiniteLength
                return Vector{Atom}(), false, 0.0 # means find nothing
            end 
            index = neighborInfo.index
            ChangeCell!(atom, index, simulator)
            cell = GetCell(grid, index, simulator)
        end
    end
end

function _append_neighbor_candidates!(
    buf::Vector{Atom},
    atom::Atom,
    neighborCellInfo::NeighborCellInfo,
    filterIndexes::Vector{Int64},
    simulator::Simulator,
)
    box = simulator.box
    pMax = simulator.parameters.pMax
    cross = neighborCellInfo.cross
    for d in 1:3
        if cross[d] != 0 && !simulator.parameters.periodic[d]
            return nothing
        end
    end
    neighborCell = GetCell(simulator.grid, neighborCellInfo.index, simulator)
    for neighborAtom in neighborCell.atoms
        if neighborAtom.index == atom.index || neighborAtom.index in filterIndexes
            continue
        end
        if ComputeVDistance(atom, neighborAtom, cross, box) > 0
            p = ComputeP!(atom, neighborAtom, cross, box)
            if p < pMax
                push!(buf, neighborAtom)
            end
        end
    end
    if neighborCell.isNonLatticeAtoms
        RefillLatticeAtoms!(neighborCell, simulator)
    end
    for neighborAtom in neighborCell.latticeAtoms
        if !neighborAtom.isAlive || neighborAtom.index in filterIndexes
            continue
        end
        if ComputeVDistance(atom, neighborAtom, cross, box) > 0
            p = ComputeP!(atom, neighborAtom, cross, box)
            if p < pMax
                push!(buf, neighborAtom)
            end
        end
    end
    return nothing
end

function _append_threaded_neighbor_candidates!(
    threadCandidates::Vector{Vector{Atom}},
    atom::Atom,
    neighborCellsInfo::Array{NeighborCellInfo, 3},
    filterIndexes::Vector{Int64},
    simulator::Simulator,
)
    grid = simulator.grid
    for n in eachindex(neighborCellsInfo)
        neighborCellInfo = neighborCellsInfo[n]
        GetCell(grid, neighborCellInfo.index, simulator)  # preload the cell to avoid race condition
    end
    @threads :static for n in eachindex(neighborCellsInfo)
        neighborCellInfo = neighborCellsInfo[n]
        buf = threadCandidates[Threads.threadid()]
        _append_neighbor_candidates!(buf, atom, neighborCellInfo, filterIndexes, simulator)
    end
    return nothing
end

function GetTargetsFromNeighbor_dynamicLoad(atom::Atom, cell::Cell, filterIndexes::Vector{Int64}, simulator::Simulator)
    targets = Vector{Atom}()
    nthreads = Threads.nthreads()
    neighborCellsInfo = GetNeighborCellsInfo!(cell, simulator.grid, simulator)
    threadCandidates = simulator.workBuffers.threadCandidates
    for tc in threadCandidates
        empty!(tc)
    end
    if nthreads == 1
        buf = threadCandidates[1]
        for neighborCellInfo in neighborCellsInfo
            _append_neighbor_candidates!(buf, atom, neighborCellInfo, filterIndexes, simulator)
        end
    else
        _append_threaded_neighbor_candidates!(threadCandidates, atom, neighborCellsInfo, filterIndexes, simulator)
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
    push!(simulator.preservedCellKeys, nearestTarget.cellIndex)

    for candidateTarget in candidateTargets
        if candidateTarget.index == nearestTarget.index
            continue
        end
        if SimultaneousCriteria(candidateTarget, nearestTarget, simulator)
            push!(targets, candidateTarget)
            push!(simulator.preservedCellKeys, candidateTarget.cellIndex)
        end
    end    
    return targets
end


function Collision_dynamicLoad!(atom_p::Atom, atoms_t::Vector{Atom}, emptyPath::Float64, simulator::Simulator)
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
    pL -= emptyPath
    N = simulator.uniformDensity
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
            simulator.θFunctions[(atom_p.type, atom_t.type)], simulator.τFunctions[(atom_p.type, atom_t.type)])
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
    for i in eachindex(E_tList)
        E_tList[i] *= λ
    end
    for (i, atom_t) in enumerate(atoms_t)
        if E_tList[i] > GetDTE(atom_t, simulator) && E_tList[i] - GetBDE(atom_t, simulator) > 0.1
            SetEnergy!(atom_t, E_tList[i] - GetBDE(atom_t, simulator))
        else
            SetEnergy!(atom_t, 0.0)
        end
    end
end 




function Cascade_dynamicLoad!(atom_p::Atom, simulator::Simulator)
    pAtoms = Atom[atom_p]
    pAtomsIndex = Int64[atom_p.index]
    filterIndexes = Int64[]
    parameters = simulator.parameters
    simulator.nCollisionEvent = 0
    simulator.nCascade += 1
    DumpInCascade_dynamicLoad(simulator)
    while true
        simulator.nCollisionEvent += 1
        targetsList = Vector{Vector{Atom}}()
        emptyPathList = Float64[]
        deleteIndexes = Int64[]
        othersTargetIndexes = Int64[]
        for (na, pAtom) in enumerate(pAtoms)
            empty!(filterIndexes)
            append!(filterIndexes, pAtomsIndex)
            append!(filterIndexes, LastTargets!(pAtom, simulator))
            append!(filterIndexes, othersTargetIndexes)
            targets, isAlive, emptyPath = ShotTarget_dynamicLoad(pAtom, filterIndexes, simulator)
            if !isAlive
                ClearLastTargets!(pAtom, simulator)
                delete_dynamicLoad!(simulator, pAtom)
                push!(deleteIndexes, na)
                continue
            end
            push!(targetsList, targets)
            push!(emptyPathList, emptyPath)
            for target in targets
                push!(othersTargetIndexes, target.index)
            end
        end
        deleteat!(pAtoms, deleteIndexes)
        empty!(pAtomsIndex)
        for pAtom in pAtoms
            push!(pAtomsIndex, pAtom.index)
        end
        nextPAtoms = Vector{Atom}()
        for (pAtom, targets, emptyPath) in zip(pAtoms, targetsList, emptyPathList)
            if length(targets) > 0
                lastTargets = LastTargets!(pAtom, simulator)
                empty!(lastTargets)
                for target in targets
                    push!(lastTargets, target.index)
                end
                Collision_dynamicLoad!(pAtom, targets, emptyPath, simulator)
                for target in targets
                    if target.energy > 0.0   
                        if target.isLatticeAtom
                            target = LeaveLatticePoint_dynamicLoad!(target, simulator)
                        end
                        #DisplaceAtom!(target, target.coordinate, simulator) # why I do this?
                        push!(nextPAtoms, target)
                        targetLastTargets = LastTargets!(target, simulator)
                        empty!(targetLastTargets)
                        push!(targetLastTargets, pAtom.index)
                    end
                end
                if pAtom.energy > parameters.stopEnergy 
                    push!(nextPAtoms, pAtom)
                else
                    ClearLastTargets!(pAtom, simulator)
                    Stop_dynamicLoad!(pAtom, simulator)
                end
            else
                push!(nextPAtoms, pAtom)
            end
        end
        empty!(simulator.preservedCellKeys)
        for k in simulator.attempedDeCellKeys
            DeprecateCell!(GetCell(simulator.grid, k, simulator), simulator)
        end
        empty!(simulator.attempedDeCellKeys)
        DumpInCascade_dynamicLoad(simulator)
        if length(nextPAtoms) > 0
            pAtoms = nextPAtoms
            sort!(pAtoms, by = a -> a.energy, rev = true)
            empty!(pAtomsIndex)
            for pAtom in pAtoms
                push!(pAtomsIndex, pAtom.index)
            end
        else
            break
        end
    end
    empty!(simulator.workBuffers.lastTargets)
end




function delete_dynamicLoad!(simulator::Simulator, atom::Atom; isDeleteVacancy::Bool = false)
    cell = GetCell(simulator.grid, atom.cellIndex, simulator)
    if !isDeleteVacancy
        deleteat!(cell.atoms, findfirst(a -> a.index == atom.index, cell.atoms))
        simulator.numberOfAtoms -= 1
    else
        deleteat!(cell.vacancies, findfirst(v -> v.index == atom.index, cell.vacancies))
        cell.latticeAtoms[atom.indexInCell].isAlive = true
        simulator.numberOfVacancies -= 1
    end
    atom.isAlive = false 
    DeprecateCell!(cell, simulator)
end


function Stop_dynamicLoad!(atom::Atom, simulator::Simulator)
    grid = simulator.grid
    cell = GetCell(grid, atom.cellIndex, simulator)
    nearestVacancyDistance_squared = Inf
    isExist = false
    nearestVacancy = nothing  
    nearestCell = nothing
    neighborCellsInfo = GetNeighborCellsInfo!(cell, grid, simulator)
    for neighborCellInfo in neighborCellsInfo
        index = neighborCellInfo.index
        cross = neighborCellInfo.cross
        neighborCell = GetCell(simulator.grid, index, simulator)
        DeprecateCell!(neighborCell, simulator)
        if simulator.parameters.vacancyRecoverDistance_squared == 0.0
            continue
        end
        for vacancy in neighborCell.vacancies
            dr2 = ComputeDistance_squared(atom.coordinate, vacancy.coordinate, cross, simulator.box)
            if dr2 < simulator.parameters.vacancyRecoverDistance_squared && dr2 < nearestVacancyDistance_squared
                nearestVacancyDistance_squared = dr2
                nearestVacancy = vacancy  # store the nearest vacancy
                nearestCell = neighborCell
                isExist = true
            end
        end
    end
    if isExist && nearestVacancy !== nothing
        if atom.type == nearestVacancy.type - length(keys(simulator.parameters.typeDict)) 
            nearestCell.latticeAtoms[nearestVacancy.indexInCell].isAlive = true
            delete_dynamicLoad!(simulator, atom)
            delete_dynamicLoad!(simulator, nearestVacancy, isDeleteVacancy = true)
        else
            SetCoordinate!(atom, nearestVacancy.coordinate)
            Pertubation_dynamicload!(atom, nearestCell.ranges, simulator)
            ChangeCell!(atom, nearestVacancy.cellIndex, simulator)
        end
    end
end


function LeaveLatticePoint_dynamicLoad!(latticeAtom::Atom, simulator::Simulator; isUpdateEnv::Bool = true)
    cell = GetCell(simulator.grid, latticeAtom.cellIndex, simulator)
    vacancy = CreateVacancy(latticeAtom, simulator)
    push!(cell.vacancies, vacancy)
    push!(simulator.vacancies, vacancy)
    simulator.numberOfVacancies += 1
    vacancy.index = simulator.maxVacancyID
    simulator.maxVacancyID += 1
    vacancy.indexInCell = latticeAtom.indexInCell
    vacancy.cellIndex = latticeAtom.cellIndex

    cell.latticeAtoms[latticeAtom.indexInCell].isAlive = false

    atom = Atom(latticeAtom.type, latticeAtom.coordinate, simulator.parameters)
    atom.isLatticeAtom = false
    atom.velocityDirection = latticeAtom.velocityDirection
    atom.energy = latticeAtom.energy
    for d in 1:3
        length = simulator.box.vectors[d,d]
        if atom.coordinate[d] < 0
            atom.coordinate[d] += length
        elseif atom.coordinate[d] >= length
            atom.coordinate[d] -= length
        end
    end
    cellIndex = WhichCell(atom.coordinate, simulator.grid)
    if latticeAtom.cellIndex != cellIndex
        DeprecateCell!(cell, simulator)
        cell = GetCell(simulator.grid,cellIndex, simulator)
    end
    push!(cell.atoms, atom)
    atom.cellIndex = cell.index
    simulator.maxAtomID += 1
    simulator.numberOfAtoms += 1
    atom.index = simulator.maxAtomID
    push!(simulator.atoms, atom)
    return atom
end

function CreateVacancy(atom::Atom, simulator::Simulator)
    vacancy = Atom(atom.type, atom.latticeCoordinate, simulator.parameters)
    vacancy.cellIndex = atom.cellIndex
    vacancy.type += length(keys(simulator.parameters.typeDict))
    vacancy.indexInCell = atom.indexInCell
    return vacancy
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
    empty!(simulator.grid.cells)
    empty!(simulator.deprecatedCellKeys)
    empty!(simulator.preservedCellKeys)
    empty!(simulator.attempedDeCellKeys)
    empty!(simulator.atoms)
    empty!(simulator.vacancies)
    simulator.maxAtomID = 0
    simulator.maxVacancyID = 1E6 
    simulator.minLatticeAtomID = 0
    simulator.numberOfAtoms = 0
    simulator.numberOfVacancies = 0
end


function Pertubation_dynamicload!(atom::Atom, ranges::Matrix{Float64}, simulator::Simulator)
    if simulator.parameters.isAmorphous 
        rng = THREAD_RNG[Threads.threadid()]
        atom.coordinate[1] = ranges[1, 1] + rand(rng) * simulator.grid.vectors[1, 1]
        atom.coordinate[2] = ranges[2, 1] + rand(rng) * simulator.grid.vectors[2, 2]
        atom.coordinate[3] = ranges[3, 1] + rand(rng) * simulator.grid.vectors[3, 3]
    else
        ah = simulator.parameters.amorphousHeight
        if atom.coordinate[3] > ah
            rng = THREAD_RNG[Threads.threadid()]
            atom.coordinate[1] = ranges[1,1] + rand(rng) * simulator.grid.vectors[1, 1]
            atom.coordinate[2] = ranges[2,1] + rand(rng) * simulator.grid.vectors[2, 2]
            base = ranges[3,1] > ah ? ranges[3,1] : ah
            latticeTop = simulator.parameters.primaryVectors[3,3] * simulator.parameters.latticeRanges[3,2]     
            top = ranges[3,2] < latticeTop ? ranges[3,2] : latticeTop
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

function DumpInCascade_dynamicLoad(simulator::Simulator)
    if simulator.parameters.isDumpInCascade
        if simulator.parameters.debugMode == false
            @dump "Cascade_$(simulator.nCascade).dump" [simulator.atoms; simulator.vacancies] ["vx", "vy", "vz", "e"]
        else
            cells = values(simulator.grid.cells)
            atoms = Vector{Atom}()
            for cell in cells
                if !(cell.index in simulator.deprecatedCellKeys)
                    append!(atoms, [atom for atom in cell.latticeAtoms if atom.isAlive])
                    append!(atoms, cell.atoms)
                end
            end
            @dump "Cascade_$(simulator.nCascade).dump" atoms ["vx", "vy", "vz", "e", "isLatticeAtom"]
        end
    end
end
