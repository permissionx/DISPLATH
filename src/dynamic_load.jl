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
            dimension, direction, t = AtomOutFaceDimension(atom, cell, simulator)
            emptyPath = t
            ni = direction == 1 ? 1 : 3
            neighborCellsInfo = GetNeighborCellsInfo!(cell, grid, simulator)
            neighborInfo = neighborCellsInfo[dimension == 1 ? ni : 2,
                                             dimension == 2 ? ni : 2,
                                             dimension == 3 ? ni : 2]
            crossFlag = neighborInfo.cross
            if crossFlag[dimension] != 0 && periodic[dimension]
                atom.coordinate[dimension] -= crossFlag[dimension] * simulator.box.vectors[dimension, dimension]
            end
            if (neighborInfo.cross[dimension] != 0 && !periodic[dimension]) || t >= simulator.parameters.infiniteLength
                return TargetCandidate[], false, 0.0 # means find nothing
            end
            index = neighborInfo.index
            ChangeCell!(atom, index, simulator)
            cell = GetCell(grid, index, simulator)
        end
    end
end

function _append_neighbor_candidates!(
    buf::Vector{TargetCandidate},
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
        if ComputeVDistance(atom, neighborAtom, cross, box, simulator) > 0
            candidate = ComputeP(atom, neighborAtom, cross, box, simulator)
            if candidate.pValue < pMax
                push!(buf, candidate)
            end
        end
    end
    if IsEmptyDynamicCell(neighborCell.index, simulator)
        return nothing
    end
    latticeCoordinates = LatticeSiteCoordinates!(neighborCell.index, simulator)
    for (indexInCell, stdAtom) in enumerate(simulator.cellStd.atoms)
        targetIndex = LatticeSiteIndex(neighborCell.index, indexInCell, simulator)
        if targetIndex in filterIndexes || HasVacancyAtIndex(neighborCell, indexInCell, simulator)
            continue
        end
        coordinate = latticeCoordinates[indexInCell]
        if ComputeVDistance(atom, coordinate, cross, box, simulator) > 0
            candidate = ComputeP(atom, targetIndex, stdAtom.type, neighborCell.index, true, indexInCell, coordinate, cross, box, simulator)
            if candidate.pValue < pMax
                push!(buf, candidate)
            end
        end
    end
    return nothing
end

function _append_threaded_neighbor_candidates!(
    threadCandidates::Vector{Vector{TargetCandidate}},
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
    targets = Vector{TargetCandidate}()
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


function Collision_dynamicLoad!(atom_p::Atom, targets::Vector{TargetCandidate}, emptyPath::Float64, simulator::Simulator)
    N_t = length(targets)
    grid = simulator.grid
    buffers = simulator.workBuffers.collisionParames
    EnsureCollisionCapacity!(buffers, N_t)
    tanφList = @view buffers.tanφList[1:N_t]
    tanψList = @view buffers.tanψList[1:N_t]
    E_tList = @view buffers.E_tList[1:N_t]
    x_pList = @view buffers.x_pList[1:N_t]
    x_tList = @view buffers.x_tList[1:N_t]
    Q_locList = @view buffers.Q_locList[1:N_t]
    firstTarget = targets[1]
    pL = firstTarget.pL
    pPoint = firstTarget.pPoint
    pL -= emptyPath
    N = simulator.uniformDensity
    atom_p_energy = AtomEnergy(atom_p, simulator)
    atom_p_mass = AtomMass(atom_p, simulator)
    Q_nl_v = Q_nl(atom_p_energy, atom_p_mass, TargetMass(firstTarget, simulator), atom_p.type, firstTarget.type,
                         pL, N, simulator.constantsByType)
    atom_p_energy -= Q_nl_v
    #if atom_p.type == 2
    #global Q_loss += Q_nl_v  # debug 
    #end
    if atom_p_energy < 0.1 && atom_p_energy + Q_nl_v >= 0.1
        atom_p_energy = 0.11
    end
    momentum = @SVector [0.0, 0.0, 0.0] 
    atom_p_velocity = AtomVelocityDirection(atom_p, simulator)
    for (i, target) in enumerate(targets)
        p = target.pValue
        #N = simulator.uniformDensity 
        tanφList[i], tanψList[i], E_tList[i], x_pList[i], x_tList[i], Q_locList[i] = CollisionParams(
            atom_p_energy, atom_p_mass, TargetMass(target, simulator), atom_p.type, target.type, p, simulator.constantsByType,
            simulator.θFunctions[(atom_p.type, target.type)], simulator.τFunctions[(atom_p.type, target.type)])
        if target.pValue != 0
            velocityDirectionTmp = -target.pVector / target.pValue * tanψList[i] + atom_p_velocity
        else
            velocityDirectionTmp = atom_p_velocity
        end   
        SetVelocityDirection!(target, velocityDirectionTmp, simulator)
        momentum += sqrt(2 * TargetMass(target, simulator) * E_tList[i]) * AtomVelocityDirection(target, simulator)
    end
    pMomentum = sqrt(2 * atom_p_mass * atom_p_energy) * atom_p_velocity - momentum
    pVelocity = pMomentum / atom_p_mass
    SetVelocityDirection!(atom_p, pVelocity, simulator)
    pEnergy =  sum(pMomentum .* pMomentum) / 2 / atom_p_mass
    sumE_t = sum(E_tList)
    sumQ_loc = sum(Q_locList) 
    ENeed = atom_p_energy - sumQ_loc # - (N_t - 1) * Q_nl_v
    λ = ENeed / (pEnergy + sumE_t)
    DisplaceAtom!(atom_p, pPoint, simulator)
    SetEnergy!(atom_p, pEnergy * λ, simulator)
    #if atom_p.type == 2
    #    @record "log/$(simulator.nCascade).csv" "$(pEnergy * λ),$(minimum([a.pValue for a in atoms_t])),$(pL),$(N_t),$(atom_p.coordinate[1]),$(atom_p.coordinate[2]),$(atom_p.coordinate[3]),$(atom_p.velocityDirection[1]),$(atom_p.velocityDirection[2]),$(atom_p.velocityDirection[3])" "e,p,pL,N_t,x,y,z,vx,vy,vz" 
    #end
    for i in eachindex(E_tList)
        E_tList[i] *= λ
    end
    for (i, target) in enumerate(targets)
        if E_tList[i] > GetDTE(target, simulator) && E_tList[i] - GetBDE(target, simulator) > 0.1
            SetEnergy!(target, E_tList[i] - GetBDE(target, simulator), simulator)
        else
            SetEnergy!(target, 0.0, simulator)
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
    ClearLatticeSiteCoordinateCaches!(simulator.workBuffers)
    DumpInCascade_dynamicLoad(simulator)
    while true
        simulator.nCollisionEvent += 1
        targetsList = Vector{Vector{TargetCandidate}}()
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
                for targetCandidate in targets
                    push!(lastTargets, targetCandidate.index)
                end
                Collision_dynamicLoad!(pAtom, targets, emptyPath, simulator)
                for targetCandidate in targets
                    if AtomEnergy(targetCandidate, simulator) > 0.0
                        if targetCandidate.isLatticeAtom
                            target = LeaveLatticePoint_dynamicLoad!(targetCandidate, simulator)
                        else
                            target = TargetAtom(targetCandidate, simulator)
                        end
                        #DisplaceAtom!(target, target.coordinate, simulator) # why I do this?
                        push!(nextPAtoms, target)
                        targetLastTargets = LastTargets!(target, simulator)
                        empty!(targetLastTargets)
                        push!(targetLastTargets, pAtom.index)
                    end
                end
                if AtomEnergy(pAtom, simulator) > parameters.stopEnergy
                    push!(nextPAtoms, pAtom)
                else
                    ClearLastTargets!(pAtom, simulator)
                    Stop_dynamicLoad!(pAtom, simulator)
                end
            else
                push!(nextPAtoms, pAtom)
            end
        end
        DumpInCascade_dynamicLoad(simulator)
        if length(nextPAtoms) > 0
            pAtoms = nextPAtoms
            sort!(pAtoms, by = a -> AtomEnergy(a, simulator), rev = true)
            empty!(pAtomsIndex)
            for pAtom in pAtoms
                push!(pAtomsIndex, pAtom.index)
            end
        else
            break
        end
    end
    empty!(simulator.workBuffers.lastTargets)
    empty!(simulator.workBuffers.atomDynamics)
    ClearLatticeSiteCoordinateCaches!(simulator.workBuffers)
    SweepTouchedCells!(simulator)
end




function delete_dynamicLoad!(simulator::Simulator, atom::Atom; isDeleteVacancy::Bool = false)
    cell = GetCell(simulator.grid, atom.cellIndex, simulator)
    if !isDeleteVacancy
        deleteat!(cell.atoms, findfirst(a -> a.index == atom.index, cell.atoms))
        simulator.numberOfAtoms -= 1
    else
        deleteat!(cell.vacancies, findfirst(v -> v.index == atom.index, cell.vacancies))
        simulator.numberOfVacancies -= 1
    end
    atom.isAlive = false
    ClearAtomDynamics!(atom, simulator)
    MarkCellIfEmpty!(cell, simulator)
end


function Stop_dynamicLoad!(atom::Atom, simulator::Simulator)
    ClearAtomDynamics!(atom, simulator)
    if simulator.parameters.vacancyRecoverDistance_squared == 0.0
        return
    end
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
            delete_dynamicLoad!(simulator, atom)
            delete_dynamicLoad!(simulator, nearestVacancy, isDeleteVacancy = true)
        else
            SetCoordinate!(atom, nearestVacancy.coordinate)
            Pertubation_dynamicload!(atom, CellRanges(nearestCell, simulator.grid), simulator)
            ChangeCell!(atom, nearestVacancy.cellIndex, simulator)
        end
    end
end

function LeaveLatticePoint_dynamicLoad!(target::TargetCandidate, simulator::Simulator; isUpdateEnv::Bool = true)
    cell = GetCell(simulator.grid, target.cellIndex, simulator)
    vacancy = CreateVacancy(target, simulator)
    push!(cell.vacancies, vacancy)
    push!(simulator.vacancies, vacancy)
    simulator.numberOfVacancies += 1
    vacancy.index = simulator.maxVacancyID
    simulator.maxVacancyID += 1

    velocityDirection = AtomVelocityDirection(target, simulator)
    energy = AtomEnergy(target, simulator)
    ClearAtomDynamics!(target.index, simulator)
    atom = Atom(target.type, target.coordinate, simulator.parameters)
    for d in 1:3
        length = simulator.box.vectors[d,d]
        if atom.coordinate[d] < 0
            atom.coordinate[d] += length
        elseif atom.coordinate[d] >= length
            atom.coordinate[d] -= length
        end
    end
    cellIndex = WhichCell(atom.coordinate, simulator.grid)
    if target.cellIndex != cellIndex
        cell = GetCell(simulator.grid, cellIndex, simulator)
    end
    push!(cell.atoms, atom)
    atom.cellIndex = cell.index
    simulator.maxAtomID += 1
    simulator.numberOfAtoms += 1
    atom.index = simulator.maxAtomID
    push!(simulator.atoms, atom)
    SetVelocityDirection!(atom, velocityDirection, simulator)
    SetEnergy!(atom, energy, simulator)
    return atom
end


function LeaveLatticePoint_dynamicLoad!(latticeAtom::Atom, simulator::Simulator; isUpdateEnv::Bool = true)
    error("Dynamic lattice atoms are represented by TargetCandidate, not stored Atom objects.")
end

function CreateVacancy(atom::Atom, simulator::Simulator)
    vacancy = Atom(atom.type, LatticeCoordinate(atom, simulator), simulator.parameters)
    vacancy.cellIndex = atom.cellIndex
    vacancy.type += length(keys(simulator.parameters.typeDict))
    return vacancy
end

function CreateVacancy(target::TargetCandidate, simulator::Simulator)
    vacancy = Atom(target.type, LatticeCoordinate(target.cellIndex, target.indexInCell, simulator), simulator.parameters)
    vacancy.cellIndex = target.cellIndex
    vacancy.type += length(keys(simulator.parameters.typeDict))
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
                velocityDirection = AtomVelocityDirection(atom, simulator)
                energy = AtomEnergy(atom, simulator)
                mass = AtomMass(atom, simulator)
                if isDebug
                    write(file, "$(atom.index) $(atom.type) \
                    $(atom.coordinate[1]) $(atom.coordinate[2]) $(atom.coordinate[3]) \
                    $(velocityDirection[1]*sqrt(2*mass*energy)) $(velocityDirection[2]*sqrt(2*mass*energy)) $(velocityDirection[3]*sqrt(2*mass*energy)) \
                    $(energy) \
                    $(atom.cellIndex[1]) $(atom.cellIndex[2]) $(atom.cellIndex[3]) \
                    $(GetDTE(atom, simulator))\n")
                else
                    write(file, "$(atom.index) $(atom.type) \
                    $(atom.coordinate[1]) $(atom.coordinate[2]) $(atom.coordinate[3]) $(energy)\n")
                end
            end
        end 
        for atom in simulator.vacancies
            if atom.isAlive
                velocityDirection = AtomVelocityDirection(atom, simulator)
                energy = AtomEnergy(atom, simulator)
                mass = AtomMass(atom, simulator)
                if isDebug
                    write(file, "$(atom.index+100000) $(atom.type) \
                    $(atom.coordinate[1]) $(atom.coordinate[2]) $(atom.coordinate[3]) \
                    $(velocityDirection[1]*sqrt(2*mass*energy)) $(velocityDirection[2]*sqrt(2*mass*energy)) $(velocityDirection[3]*sqrt(2*mass*energy)) \
                    $(energy) \
                    $(atom.cellIndex[1]) $(atom.cellIndex[2]) $(atom.cellIndex[3]) \
                    $(GetDTE(atom, simulator))\n")
                else
                    write(file, "$(atom.index+100000) $(atom.type) \
                    $(atom.coordinate[1]) $(atom.coordinate[2]) $(atom.coordinate[3]) $(energy)\n")
                end
            end
        end
    end
end



function Restore_dynamicLoad!(simulator::Simulator)
    empty!(simulator.grid.cells)
    empty!(simulator.freeCells)
    empty!(simulator.touchedCells)
    empty!(simulator.atoms)
    empty!(simulator.vacancies)
    simulator.maxAtomID = 0
    simulator.maxVacancyID = 1E6 
    simulator.minLatticeAtomID = 0
    simulator.numberOfAtoms = 0
    simulator.numberOfVacancies = 0
    ClearBuffers!(simulator.workBuffers)
end


function Pertubation_dynamicload!(atom::Atom, ranges::AbstractMatrix{<:Real}, simulator::Simulator)
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
                if !IsEmptyDynamicCell(cell.index, simulator)
                    for (indexInCell, stdAtom) in enumerate(simulator.cellStd.atoms)
                        HasVacancyAtIndex(cell, indexInCell, simulator) && continue
                        atom = Atom(stdAtom.type, LatticeSiteCoordinate(cell.index, indexInCell, simulator), simulator.parameters)
                        atom.index = LatticeSiteIndex(cell.index, indexInCell, simulator)
                        atom.cellIndex = cell.index
                        push!(atoms, atom)
                    end
                end
                append!(atoms, cell.atoms)
            end
            @dump "Cascade_$(simulator.nCascade).dump" atoms ["vx", "vy", "vz", "e", "isLatticeAtom"]
        end
    end
end
