using StaticArrays

function Box(Vectors::Matrix{Float64})
    log_info("Box created: $(round(Vectors[1,1]; digits=2)) × $(round(Vectors[2,2]; digits=2)) × $(round(Vectors[3,3]; digits=2)) Å")
    return Box(Vectors, inv(Vectors'), true)
end 

function CreateBoxByPrimaryVectors(primaryVectors::Matrix{Float64}, sizes::Vector{Int64})
    vectors = primaryVectors .* sizes
    return Box(vectors)
end 


function Atom(type::Int64, coordinate::AbstractVector{<:Real}, parameters::Parameters)
    index = 0
    isAlive = true
    cellIndex = (0,0,0)
    if IS_DYNAMIC_LOAD
        return Atom(index, isAlive, type,
                    SVector{3,Float64}(coordinate[1], coordinate[2], coordinate[3]), cellIndex)
    end
    coordinateVector = Float64[coordinate[1], coordinate[2], coordinate[3]]

    velocityDirection = ZERO_VELOCITY_DIRECTION
    energy = 0.0
    radius, mass, Z, dte, bde, _, _ = TypeToProperties(type, parameters.typeDict)
    emptyPath = 0.0
    pValue = 0.0
    pPoint = SVector{3,Float64}(0.0, 0.0, 0.0)
    pVector = SVector{3,Float64}(0.0, 0.0, 0.0)
    pL = 0.0
    pAtomIndex = -1
    pDirection = Float64[0.0, 0.0, 0.0]
    lastTargets = Vector{Int64}()
    latticePointIndex = -1
    frequency = 0.0
    frequencies = Vector{Float64}()
    finalLatticePointIndexs = Vector{Int64}()
    eventIndex = -1
    isNewlyLoaded = false
    latticeCoordinate = SVector{3,Float64}(coordinateVector[1], coordinateVector[2], coordinateVector[3])
    indexInCell = 0
    return Atom(index, isAlive, type, coordinateVector, cellIndex,
                radius, mass, velocityDirection, energy, Z,
                dte, bde, emptyPath,
                pValue, pPoint, pVector, pL, pAtomIndex, pDirection,
                lastTargets, latticePointIndex,
                frequency, frequencies, finalLatticePointIndexs, eventIndex,
                isNewlyLoaded, latticeCoordinate, indexInCell)
end


function TypeToProperties(type::Int64, typeDict::Dict{Int64, Element})
    if haskey(typeDict, type)
        element = typeDict[type]
        return element.radius, element.mass, element.Z, element.dte, element.bde, element.alpha, element.beta 
    else
        error("Unknown atom type: $type")
    end 
end 

const PENDING_ATOM_DYNAMICS = IdDict{Atom, AtomDynamics}()

@inline function BaseType(type::Int64, simulator::Simulator)
    ntypes = length(simulator.parameters.typeDict)
    return type > ntypes ? type - ntypes : type
end

@inline function TypeElement(type::Int64, simulator::Simulator)
    return simulator.parameters.typeDict[BaseType(type, simulator)]
end

@inline function AtomElement(atom::Atom, simulator::Simulator)
    return TypeElement(atom.type, simulator)
end

@inline AtomMass(atom::Atom, simulator::Simulator) = AtomElement(atom, simulator).mass
@inline AtomDTE(atom::Atom, simulator::Simulator) = AtomElement(atom, simulator).dte
@inline AtomBDE(atom::Atom, simulator::Simulator) = AtomElement(atom, simulator).bde
@inline TargetMass(target::TargetCandidate, simulator::Simulator) = TypeElement(target.type, simulator).mass
@inline TargetDTE(target::TargetCandidate, simulator::Simulator) = TypeElement(target.type, simulator).dte
@inline TargetBDE(target::TargetCandidate, simulator::Simulator) = TypeElement(target.type, simulator).bde

function _pending_dynamics(atom::Atom)
    return get(PENDING_ATOM_DYNAMICS, atom, ZERO_ATOM_DYNAMICS)
end

function _transfer_pending_dynamics!(atom::Atom, simulator::Simulator)
    if !IS_DYNAMIC_LOAD
        return nothing
    end
    dynamics = get(PENDING_ATOM_DYNAMICS, atom, nothing)
    if dynamics !== nothing
        if dynamics.energy > 0.0
            simulator.workBuffers.atomDynamics[atom.index] = dynamics
        end
        delete!(PENDING_ATOM_DYNAMICS, atom)
    end
    return nothing
end

function AtomEnergy(atom::Atom, simulator::Simulator)
    if !IS_DYNAMIC_LOAD
        return atom.energy
    end
    return AtomEnergy(atom.index, simulator)
end

function AtomEnergy(index::Int64, simulator::Simulator)
    return get(simulator.workBuffers.atomDynamics, index, ZERO_ATOM_DYNAMICS).energy
end

function AtomEnergy(target::TargetCandidate, simulator::Simulator)
    return AtomEnergy(target.index, simulator)
end

function AtomVelocityDirection(atom::Atom, simulator::Simulator)
    if !IS_DYNAMIC_LOAD
        return atom.velocityDirection
    end
    return AtomVelocityDirection(atom.index, simulator)
end

function AtomVelocityDirection(index::Int64, simulator::Simulator)
    return get(simulator.workBuffers.atomDynamics, index, ZERO_ATOM_DYNAMICS).velocityDirection
end

function AtomVelocityDirection(target::TargetCandidate, simulator::Simulator)
    return AtomVelocityDirection(target.index, simulator)
end

function CreateGrid(box::Box, inputVectors::Matrix{Float64})
    if !box.isOrthogonal
        error("The box is not orthogonal, please use the orthogonal box.")
    end
    sizes = Vector{Int64}(undef, 3)
    vectors = zeros(Float64, 3, 3)

    if ! IS_DYNAMIC_LOAD
        for d in 1:3
            sizes[d] = Int64(floor(box.vectors[d,d] / inputVectors[d,d]))
            vectors[d,d] = box.vectors[d,d] / sizes[d]
        end
    else
        for d in 1:3
            sizes[d] = Int64(round(box.vectors[d,d] / inputVectors[d,d]))
            vectors[d,d] = inputVectors[d,d]
        end
    end
    for d in 1:3
        if sizes[d] < 3
            error("The box size in dimension $d is too small,  use a larger box! (At least 3 cells in each dimension)")
        end
    end
    log_info("Cell grid: $(sizes[1]) × $(sizes[2]) × $(sizes[3]) = $(sizes[1]*sizes[2]*sizes[3]) cells")
    log_info("Cell size: $(round(vectors[1,1]; digits=2)) × $(round(vectors[2,2]; digits=2)) × $(round(vectors[3,3]; digits=2)) Å")
    if ! IS_DYNAMIC_LOAD
        cells = Array{Cell, 3}(undef, sizes[1], sizes[2], sizes[3])
        @showprogress desc="Creating cells: " for x in 1:sizes[1]
            for y in 1:sizes[2]
                for z in 1:sizes[3]
                    cells[x, y, z] = CreateCell((x, y, z), vectors)
                end
            end    
        end
        cellVolume = vectors[1,1] * vectors[2,2] * vectors[3,3]
        grid = Grid(cells, vectors, sizes, cellVolume) 
        @showprogress desc="Pushing cell neighbors: " for cell in grid.cells
            SetNeighborCellsInfo!(cell, grid)
        end
    else
        cells = Dict{Int64, Cell}()
        cellVolume = vectors[1,1] * vectors[2,2] * vectors[3,3]
        grid = Grid(cells, vectors, sizes, cellVolume)
    end
    log_success("Cell grid created")
    log_separator()
    return grid
end


function Simulator(box::Box, atoms::Vector{Atom}, inputGridVectors::Matrix{Float64}, parameters::Parameters)
    # this is the last entrace for simulator initilization 
    log_section("Initializing Simulator")
    simulator = Simulator(box, inputGridVectors, parameters) # object creation
    if IS_DYNAMIC_LOAD
        PN = Vector{Int64}(undef, 3)
        for d in 1:3 
            num = inputGridVectors[d,d] / parameters.primaryVectors[d, d]
            if !(round(num) ≈ num)
                error("InputGridVector must be integer multiple of primaryVector!")
            end
            PN[d] = Int64(round(inputGridVectors[d,d] / parameters.primaryVectors[d, d]))
        end
        InitCellStd!(simulator, PN)
    else
        LoadAtoms!(simulator, atoms)
    end
    log_success("Simulator initialized.")
    return simulator 
end

function LoadAtoms!(simulator::Simulator, atoms::Vector{Atom})
    for atom in atoms
        push!(simulator, atom) 
        latticePoint = LatticePoint(atom)
        push!(simulator, latticePoint)
    end
    for cell in simulator.grid.cells
        cell.atomicDensity = length(cell.atoms) / simulator.grid.cellVolume
    end 
    for atom in simulator.atoms
        Pertubation!(atom, simulator)
    end
    InitLatticePointEnvronment(simulator)
    log_success("$(simulator.numberOfAtoms) atoms loaded.")
end

function CreateAtomsByPrimaryVectors(parameters::Parameters)
    primaryVectors = parameters.primaryVectors
    latticeRanges = parameters.latticeRanges
    basis = parameters.basis
    basisTypes = parameters.basisTypes
    atomNumber = (latticeRanges[1,2] - latticeRanges[1,1]) * (latticeRanges[2,2] - latticeRanges[2,1]) * (latticeRanges[3,2] - latticeRanges[3,1]) * length(basisTypes)
    atoms = Vector{Atom}(undef, atomNumber)
    n = 1
    @showprogress desc="Creating atoms ($(atomNumber)): " for x in latticeRanges[1,1]:latticeRanges[1,2]-1
        for y in latticeRanges[2,1]:latticeRanges[2,2]-1    
            for z in latticeRanges[3,1]:latticeRanges[3,2]-1
                for i in eachindex(basisTypes)
                    reducedCoordinate = Float64[x,y,z] + basis[i, :]
                    coordinate = primaryVectors' * reducedCoordinate
                    atoms[n] = Atom(basisTypes[i], coordinate, parameters)
                    n += 1
                end
            end
        end
    end
    return atoms
end


function Parameters(pMax::Float64, vacancyRecoverDistance::Float64, typeDict::Dict{Int64, Element}; kwargs...)
    # non lattice info
    primaryVectors = [1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0]
    latticeRanges = [0 1; 0 1; 0 1]
    basis = [0.0 0.0 0.0]
    basisTypes = [1]  
    parameters = Parameters(primaryVectors, latticeRanges, basisTypes, basis, pMax, vacancyRecoverDistance, typeDict; kwargs...)
    return parameters
end


function LoadAtomsAndBoxFromDataFile(fileName::String; replicate::Vector{Int64} = [1,1,1])
    xlo, xhi, ylo, yhi, zlo, zhi, types, xs, ys, zs = ReadDate(fileName, replicate)
    box = Box([xhi-xlo 0.0 0.0; 0.0 yhi-ylo 0.0; 0.0 0.0 zhi-zlo])
    atoms = Vector{Atom}(undef, length(types))
    for (n,(type, x, y, z)) in enumerate(zip(types, xs, ys, zs))
        atoms[n] = Atom(type, [x, y, z], parameters)
    end
    return box, atoms
end




function LatticePoint(atom::Atom)
    environment = Vector{Int64}() 
    return LatticePoint(copy(atom.index), copy(atom.type), 
                        copy(atom.coordinate), atom.cellIndex, environment,
                        atom.index)
end


@inline function _cell_index_1d(x::Float64, d::Int64, grid::Grid)
    cellIndex = Int64(floor(x / grid.vectors[d,d])) + 1
    if cellIndex < 1
        return 1
    elseif cellIndex > grid.sizes[d]
        return grid.sizes[d]
    end
    return cellIndex
end

function WhichCell(coordinate::AbstractVector{<:Real}, grid::Grid)
    return (
        _cell_index_1d(coordinate[1], 1, grid),
        _cell_index_1d(coordinate[2], 2, grid),
        _cell_index_1d(coordinate[3], 3, grid),
    )
end

@inline IsLatticeAtom(atom::Atom) = atom.index < 0

function TargetAtom(candidate::TargetCandidate, simulator::Simulator)
    candidate.isLatticeAtom && error("Lattice target $(candidate.index) has no stored Atom")
    cell = GetCell(simulator.grid, candidate.cellIndex, simulator)
    idx = findfirst(atom -> atom.index == candidate.index, cell.atoms)
    idx === nothing && error("Target atom $(candidate.index) is not in cell $(candidate.cellIndex)")
    return cell.atoms[idx]
end

function IndexInCell(atom::Atom, cell::Cell)
    return IndexInCellByCoordinate(atom, cell, nothing)
end

function IndexInCellByCoordinate(atom::Atom, cell::Cell, simulator)
    simulator === nothing && error("Simulator is required to infer indexInCell")
    grid = simulator.grid
    for (idx, stdAtom) in enumerate(simulator.cellStd.atoms)
        if atom.type == stdAtom.type || BaseType(atom.type, simulator) == stdAtom.type
            x = stdAtom.coordinate[1] + CellLower(cell.index, 1, grid)
            y = stdAtom.coordinate[2] + CellLower(cell.index, 2, grid)
            z = stdAtom.coordinate[3] + CellLower(cell.index, 3, grid)
            if isapprox(atom.coordinate[1], x; atol=1e-8) &&
               isapprox(atom.coordinate[2], y; atol=1e-8) &&
               isapprox(atom.coordinate[3], z; atol=1e-8)
                return idx
            end
        end
    end
    error("Could not infer indexInCell for atom $(atom.index) in cell $(cell.index)")
end

function LatticeCoordinate(atom::Atom, simulator::Simulator)
    cell = GetCell(simulator.grid, atom.cellIndex, simulator)
    indexInCell = IndexInCellByCoordinate(atom, cell, simulator)
    return LatticeCoordinate(cell.index, indexInCell, simulator)
end

function LatticeCoordinate(cellIndex::Tuple{Int64, Int64, Int64}, indexInCell::Int64, simulator::Simulator)
    stdAtom = simulator.cellStd.atoms[indexInCell]
    return SVector{3,Float64}(
        stdAtom.coordinate[1] + CellLower(cellIndex, 1, simulator.grid),
        stdAtom.coordinate[2] + CellLower(cellIndex, 2, simulator.grid),
        stdAtom.coordinate[3] + CellLower(cellIndex, 3, simulator.grid),
    )
end

function LatticeSiteIndex(cellIndex::Tuple{Int64, Int64, Int64}, indexInCell::Int64, simulator::Simulator)
    sizes = simulator.grid.sizes
    linearIndex = ((cellIndex[1] - 1) * sizes[2] + (cellIndex[2] - 1)) * sizes[3] + cellIndex[3]
    return -((linearIndex - 1) * simulator.cellLatticeAtomNumber + indexInCell)
end

@inline function _splitmix64(x::UInt64)
    x += 0x9e3779b97f4a7c15
    x = (x ⊻ (x >> 30)) * 0xbf58476d1ce4e5b9
    x = (x ⊻ (x >> 27)) * 0x94d049bb133111eb
    return x ⊻ (x >> 31)
end

@inline function _unit_random(seed::UInt64, stream::Unsigned)
    bits = _splitmix64(seed + UInt64(stream))
    return Float64(bits >> 11) * 0x1.0p-53
end

function _normal_random(seed::UInt64, stream::Unsigned)
    u1 = max(_unit_random(seed, stream), eps(Float64))
    u2 = _unit_random(seed, UInt64(stream) + 0x9e3779b97f4a7c15)
    return sqrt(-2.0 * log(u1)) * cos(2π * u2)
end

function _lattice_site_seed(cellIndex::Tuple{Int64, Int64, Int64}, indexInCell::Int64, simulator::Simulator)
    seed = UInt64(simulator.nCascade + 1)
    seed ⊻= UInt64(cellIndex[1]) * 0x9e3779b97f4a7c15
    seed ⊻= UInt64(cellIndex[2]) * 0xbf58476d1ce4e5b9
    seed ⊻= UInt64(cellIndex[3]) * 0x94d049bb133111eb
    seed ⊻= UInt64(indexInCell) * 0xd6e8feb86659fd93
    return seed
end

function LatticeSiteCoordinate(cellIndex::Tuple{Int64, Int64, Int64}, indexInCell::Int64, simulator::Simulator)
    coordinate = LatticeCoordinate(cellIndex, indexInCell, simulator)
    stdAtom = simulator.cellStd.atoms[indexInCell]
    seed = _lattice_site_seed(cellIndex, indexInCell, simulator)
    grid = simulator.grid
    lo1 = CellLower(cellIndex, 1, grid)
    lo2 = CellLower(cellIndex, 2, grid)
    lo3 = CellLower(cellIndex, 3, grid)
    if simulator.parameters.isAmorphous
        return SVector{3,Float64}(
            lo1 + _unit_random(seed, 0x01) * grid.vectors[1, 1],
            lo2 + _unit_random(seed, 0x02) * grid.vectors[2, 2],
            lo3 + _unit_random(seed, 0x03) * grid.vectors[3, 3],
        )
    elseif coordinate[3] > simulator.parameters.amorphousHeight
        hi3 = CellUpper(cellIndex, 3, grid)
        ah = simulator.parameters.amorphousHeight
        base = lo3 > ah ? lo3 : ah
        latticeTop = simulator.parameters.primaryVectors[3,3] * simulator.parameters.latticeRanges[3,2]
        top = hi3 < latticeTop ? hi3 : latticeTop
        height = max(top - base, 0.0)
        return SVector{3,Float64}(
            lo1 + _unit_random(seed, 0x01) * grid.vectors[1, 1],
            lo2 + _unit_random(seed, 0x02) * grid.vectors[2, 2],
            base + _unit_random(seed, 0x03) * height,
        )
    elseif simulator.parameters.temperature > 0.0
        sigma = simulator.constantsByType.sigma[stdAtom.type]
        return SVector{3,Float64}(
            coordinate[1] + _normal_random(seed, 0x11) * sigma,
            coordinate[2] + _normal_random(seed, 0x22) * sigma,
            coordinate[3] + _normal_random(seed, 0x33) * sigma,
        )
    end
    return coordinate
end

function LatticeSiteCoordinates!(cellIndex::Tuple{Int64, Int64, Int64}, simulator::Simulator)
    cache = simulator.workBuffers.latticeSiteCoordinates[Threads.threadid()]
    coords = get(cache, cellIndex, nothing)
    if coords === nothing
        coords = Vector{SVector{3,Float64}}(undef, simulator.cellLatticeAtomNumber)
        for indexInCell in eachindex(coords)
            coords[indexInCell] = LatticeSiteCoordinate(cellIndex, indexInCell, simulator)
        end
        cache[cellIndex] = coords
    end
    return coords
end

@inline function HasVacancyAtIndex(cell::Cell, indexInCell::Int64, simulator::Simulator)
    return (cell.vacancyMask >> (indexInCell - 1)) & UInt128(1) != 0
end


function push!(simulator::Simulator, atom::Atom)
    atom.index = simulator.maxAtomID + 1
    simulator.maxAtomID += 1
    push!(simulator.atoms, atom)
    simulator.numberOfAtoms += 1
    cellIndex = WhichCell(atom.coordinate, simulator.grid)
    atom.cellIndex = cellIndex
    if IS_DYNAMIC_LOAD
        push!(GetCell(simulator.grid, cellIndex, simulator).atoms, atom)
    else
        push!(GetCell(simulator.grid, cellIndex).atoms, atom)
    end
    _transfer_pending_dynamics!(atom, simulator)
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
    vacancy = Atom(latticePoint.type, latticePoint.coordinate, simulator.parameters)
    vacancy.index = latticePoint.index
    push!(simulator.vacancies, vacancy)
    push!(GetCell(simulator.grid, latticePoint.cellIndex).vacancies, vacancy)
    if isUpdateEnv && simulator.parameters.isKMC
        DeleteAtomEvents!(simulator, atom)
        UpdateEvents!(Set(latticePoint.environment), simulator)
    end
end


function delete!(cell::Cell, atom::Atom, simulator::Simulator)
    if !atom.isAlive
        #@show simulator.nCascade, simulator.nCollisionEvent
        error("Atom $(atom.index) is not alive when deleting")
    end
    #@show atom.index, atom.cellIndex, simulator.nCascade, simulator.nCollisionEvent, atom.coordinate
    deleteat!(cell.atoms, findfirst(a -> a.index == atom.index, cell.atoms))
    atom.cellIndex = (-1, -1, -1)
    if !IS_DYNAMIC_LOAD
        cell.atomicDensity = length(cell.atoms) / simulator.grid.cellVolume  
    end
end


@inline function _wrapped_position_component(x::Float64, d::Int64, simulator::Simulator)
    if x < 0
        if simulator.parameters.periodic[d] == false
            return 0.01
        end
        return x + simulator.box.vectors[d,d]
    elseif x >= simulator.box.vectors[d,d]
        if simulator.parameters.periodic[d] == false
            return simulator.box.vectors[d,d] - 0.01
        end
        return x - simulator.box.vectors[d,d]
    end
    return x
end

function DisplaceAtom!(atom::Atom, newPosition::Union{Vector{Float64}, SVector{3, Float64}}, simulator::Simulator)
    if IS_DYNAMIC_LOAD
        atom.coordinate = SVector{3,Float64}(
            _wrapped_position_component(newPosition[1], 1, simulator),
            _wrapped_position_component(newPosition[2], 2, simulator),
            _wrapped_position_component(newPosition[3], 3, simulator),
        )
    else
        atom.coordinate[1] = _wrapped_position_component(newPosition[1], 1, simulator)
        atom.coordinate[2] = _wrapped_position_component(newPosition[2], 2, simulator)
        atom.coordinate[3] = _wrapped_position_component(newPosition[3], 3, simulator)
    end
    cellIndex = WhichCell(atom.coordinate, simulator.grid)
    if cellIndex != atom.cellIndex
        ChangeCell!(atom, cellIndex, simulator)
    end
end


function ComputeDistance_squared(coordinate1::AbstractVector{<:Real}, coordinate2::AbstractVector{<:Real}, crossFlag::NTuple{3, Int8}, box::Box)
    dv = VectorDifference(coordinate1, coordinate2, crossFlag, box)
    distance_squared = dv[1]* dv[1] + dv[2]*dv[2] + dv[3]  * dv[3]
    return distance_squared
end


function ComputeDistance(coordinate1::AbstractVector{<:Real}, coordinate2::AbstractVector{<:Real}, crossFlag::NTuple{3, Int8}, box::Box)
    return sqrt(ComputeDistance_squared(coordinate1, coordinate2, crossFlag, box))
end


function ComputeVDistance(atom_p::Atom, atom_t::Atom, crossFlag::NTuple{3, Int8}, box::Box, simulator::Simulator)
    # v for atom_p
    dv = VectorDifference(atom_p.coordinate, atom_t.coordinate, crossFlag, box)
    return dot(dv, AtomVelocityDirection(atom_p, simulator))
end

function ComputeVDistance(atom_p::Atom, atom_t::Atom, crossFlag::NTuple{3, Int8}, box::Box)
    dv = VectorDifference(atom_p.coordinate, atom_t.coordinate, crossFlag, box)
    return dot(dv, atom_p.velocityDirection)
end

function ComputeVDistance(atom_p::Atom, targetCoordinate::SVector{3,Float64}, crossFlag::NTuple{3, Int8}, box::Box, simulator::Simulator)
    dv = VectorDifference(atom_p.coordinate, targetCoordinate, crossFlag, box)
    return dot(dv, AtomVelocityDirection(atom_p, simulator))
end


function VectorDifference(v1::AbstractVector{<:Real}, v2::AbstractVector{<:Real}, crossFlag::NTuple{3, Int8}, box::Box)
    if crossFlag == (Int8(0), Int8(0), Int8(0))
        return SVector{3,Float64}(
            v2[1] - v1[1],
            v2[2] - v1[2],
            v2[3] - v1[3],
        )
    end 
    return SVector{3,Float64}(
        v2[1] - v1[1] + crossFlag[1] * box.vectors[1,1],
        v2[2] - v1[2] + crossFlag[2] * box.vectors[2,2],
        v2[3] - v1[3] + crossFlag[3] * box.vectors[3,3]
    )
end

function ComputeP!(atom_p::Atom, atom_t::Atom, crossFlag::NTuple{3, Int8}, box::Box)
    dv = VectorDifference(atom_p.coordinate, atom_t.coordinate, crossFlag, box)
    t = dot(dv, atom_p.velocityDirection)
    atom_t.pL = t
    if 1 in crossFlag || -1 in crossFlag
        pPoint_calc = Vector{Float64}(atom_p.coordinate + t * atom_p.velocityDirection)
        for d in 1:3
            if crossFlag[d] != 0
                pPoint_calc[d] -= crossFlag[d] * box.vectors[d,d]
            end
        end
    else
        pPoint_calc = atom_p.coordinate + t * atom_p.velocityDirection
    end
    atom_t.pPoint = SVector{3,Float64}(pPoint_calc[1], pPoint_calc[2], pPoint_calc[3])
    pVector_calc = atom_t.pPoint - atom_t.coordinate
    atom_t.pVector = SVector{3,Float64}(pVector_calc[1], pVector_calc[2], pVector_calc[3])
    p = norm(atom_t.pVector)
    atom_t.pValue = p
    return p
end


function ComputeP(atom_p::Atom, atom_t::Atom, crossFlag::NTuple{3, Int8}, box::Box, simulator::Simulator)
    coordinate = SVector{3,Float64}(atom_t.coordinate[1], atom_t.coordinate[2], atom_t.coordinate[3])
    return ComputeP(atom_p, atom_t.index, atom_t.type, atom_t.cellIndex, IsLatticeAtom(atom_t), 0, coordinate, crossFlag, box, simulator)
end

# Hoisted-velocity clones of ComputeVDistance/ComputeP for the dynamic-load
# candidate search: identical bodies and function boundaries (bit-identical
# results, verified), but the projectile state comes in as arguments so the
# hot loop performs no atomDynamics lookups per candidate.
function ComputeVDistanceHoisted(pCoordinate, pVelocity::SVector{3,Float64}, targetCoordinate, crossFlag::NTuple{3, Int8}, box::Box)
    dv = VectorDifference(pCoordinate, targetCoordinate, crossFlag, box)
    return dot(dv, pVelocity)
end

function ComputePHoisted(
    pCoordinate::AbstractVector{<:Real},
    pVelocity::SVector{3,Float64},
    targetIndex::Int64,
    targetType::Int64,
    targetCellIndex::Tuple{Int64, Int64, Int64},
    isLatticeAtom::Bool,
    indexInCell::Int64,
    targetCoordinate::SVector{3,Float64},
    crossFlag::NTuple{3, Int8},
    box::Box,
)
    dv = VectorDifference(pCoordinate, targetCoordinate, crossFlag, box)
    t = dot(dv, pVelocity)
    pPoint_calc = SVector{3,Float64}(
        pCoordinate[1] + t * pVelocity[1],
        pCoordinate[2] + t * pVelocity[2],
        pCoordinate[3] + t * pVelocity[3],
    )
    if crossFlag != (Int8(0), Int8(0), Int8(0))
        pPoint_calc = SVector{3,Float64}(
            pPoint_calc[1] - crossFlag[1] * box.vectors[1,1],
            pPoint_calc[2] - crossFlag[2] * box.vectors[2,2],
            pPoint_calc[3] - crossFlag[3] * box.vectors[3,3],
        )
    end
    pVector = SVector{3,Float64}(
        pPoint_calc[1] - targetCoordinate[1],
        pPoint_calc[2] - targetCoordinate[2],
        pPoint_calc[3] - targetCoordinate[3],
    )
    p = norm(pVector)
    return TargetCandidate(targetIndex, targetType, targetCellIndex, isLatticeAtom, indexInCell, targetCoordinate, p, pPoint_calc, pVector, t)
end

function ComputeP(
    atom_p::Atom,
    targetIndex::Int64,
    targetType::Int64,
    targetCellIndex::Tuple{Int64, Int64, Int64},
    isLatticeAtom::Bool,
    indexInCell::Int64,
    targetCoordinate::SVector{3,Float64},
    crossFlag::NTuple{3, Int8},
    box::Box,
    simulator::Simulator,
)
    dv = VectorDifference(atom_p.coordinate, targetCoordinate, crossFlag, box)
    velocityDirection = AtomVelocityDirection(atom_p, simulator)
    t = dot(dv, velocityDirection)
    pPoint_calc = SVector{3,Float64}(
        atom_p.coordinate[1] + t * velocityDirection[1],
        atom_p.coordinate[2] + t * velocityDirection[2],
        atom_p.coordinate[3] + t * velocityDirection[3],
    )
    if crossFlag != (Int8(0), Int8(0), Int8(0))
        pPoint_calc = SVector{3,Float64}(
            pPoint_calc[1] - crossFlag[1] * box.vectors[1,1],
            pPoint_calc[2] - crossFlag[2] * box.vectors[2,2],
            pPoint_calc[3] - crossFlag[3] * box.vectors[3,3],
        )
    end
    pVector = SVector{3,Float64}(
        pPoint_calc[1] - targetCoordinate[1],
        pPoint_calc[2] - targetCoordinate[2],
        pPoint_calc[3] - targetCoordinate[3],
    )
    p = norm(pVector)
    # need to check periodic condition
    return TargetCandidate(targetIndex, targetType, targetCellIndex, isLatticeAtom, indexInCell, targetCoordinate, p, pPoint_calc, pVector, t)
end


function SimultaneousCriteria(candidateTarget::TargetCandidate, nearestTarget::TargetCandidate, simulator::Simulator)
    deltaPL = candidateTarget.pL - nearestTarget.pL
    if deltaPL > simulator.constantsByType.qMax[(candidateTarget.type, nearestTarget.type)]
        return false
    elseif nearestTarget.pValue * nearestTarget.pValue + deltaPL * deltaPL > simulator.parameters.pMax_squared 
        return false
    elseif candidateTarget.pValue * candidateTarget.pValue + deltaPL * deltaPL > simulator.parameters.pMax_squared 
        return false
    end
    return true
end

function SimultaneousCriteria(candidateTarget::Atom, nearestTarget::Atom, simulator::Simulator)
    deltaPL = candidateTarget.pL - nearestTarget.pL
    if deltaPL > simulator.constantsByType.qMax[(candidateTarget.type, nearestTarget.type)]
        return false
    elseif nearestTarget.pValue * nearestTarget.pValue + deltaPL * deltaPL > simulator.parameters.pMax_squared
        return false
    elseif candidateTarget.pValue * candidateTarget.pValue + deltaPL * deltaPL > simulator.parameters.pMax_squared
        return false
    end
    return true
end




function SetVelocityDirection!(atom::Atom, velocity::SVector{3,Float64}, simulator::Simulator)
    if !IS_DYNAMIC_LOAD
        return SetVelocityDirection!(atom, velocity)
    end
    SetVelocityDirection!(atom.index, velocity, simulator)
end

function SetVelocityDirection!(index::Int64, velocity::SVector{3,Float64}, simulator::Simulator)
    n = norm(velocity)
    dynamics = AtomDynamics!(index, simulator)
    nextVelocity = ZERO_VELOCITY_DIRECTION
    if isnan(n) || n == Inf || n == 0.0
        nextVelocity = ZERO_VELOCITY_DIRECTION
    else
        normalized_velocity = velocity / n
        nextVelocity = SVector{3,Float64}(normalized_velocity[1], normalized_velocity[2], normalized_velocity[3])
    end
    simulator.workBuffers.atomDynamics[index] = AtomDynamics(nextVelocity, dynamics.energy)
    return nextVelocity
end

function SetVelocityDirection!(target::TargetCandidate, velocity::SVector{3,Float64}, simulator::Simulator)
    SetVelocityDirection!(target.index, velocity, simulator)
end

function SetVelocityDirection!(atom::Atom, velocity::Vector{Float64}, simulator::Simulator)
    SetVelocityDirection!(atom, SVector{3,Float64}(velocity[1], velocity[2], velocity[3]), simulator)
end

function SetVelocityDirection!(target::TargetCandidate, velocity::Vector{Float64}, simulator::Simulator)
    SetVelocityDirection!(target, SVector{3,Float64}(velocity[1], velocity[2], velocity[3]), simulator)
end

function SetEnergy!(atom::Atom, energy::Float64, simulator::Simulator)
    if !IS_DYNAMIC_LOAD
        return SetEnergy!(atom, energy)
    end
    SetEnergy!(atom.index, energy, simulator)
end

function SetEnergy!(target::TargetCandidate, energy::Float64, simulator::Simulator)
    SetEnergy!(target.index, energy, simulator)
end

function SetEnergy!(index::Int64, energy::Float64, simulator::Simulator)
    nextEnergy = energy < 0.0 ? 0.0 : energy
    if nextEnergy == 0.0
        ClearAtomDynamics!(index, simulator)
        return nothing
    end
    dynamics = AtomDynamics!(index, simulator)
    simulator.workBuffers.atomDynamics[index] = AtomDynamics(dynamics.velocityDirection, nextEnergy)
    return nothing
end

function SetVelocityDirection!(atom::Atom, velocity::SVector{3,Float64})
    if !IS_DYNAMIC_LOAD
        n = norm(velocity)
        if isnan(n) || n == Inf || n == 0.0
            atom.velocityDirection = ZERO_VELOCITY_DIRECTION
        else
            normalized_velocity = velocity / n
            atom.velocityDirection = SVector{3,Float64}(normalized_velocity[1], normalized_velocity[2], normalized_velocity[3])
        end
        return nothing
    end
    dynamics = _pending_dynamics(atom)
    n = norm(velocity)
    nextVelocity = ZERO_VELOCITY_DIRECTION
    if isnan(n) || n == Inf || n == 0.0
        nextVelocity = ZERO_VELOCITY_DIRECTION
    else
        normalized_velocity = velocity / n
        nextVelocity = SVector{3,Float64}(normalized_velocity[1], normalized_velocity[2], normalized_velocity[3])
    end
    PENDING_ATOM_DYNAMICS[atom] = AtomDynamics(nextVelocity, dynamics.energy)
end

function SetVelocityDirection!(atom::Atom, velocity::Vector{Float64})
    SetVelocityDirection!(atom, SVector{3,Float64}(velocity[1], velocity[2], velocity[3]))
end

function SetEnergy!(atom::Atom, energy::Float64)
    if !IS_DYNAMIC_LOAD
        atom.energy = energy < 0.0 ? 0.0 : energy
        return nothing
    end
    nextEnergy = energy < 0.0 ? 0.0 : energy
    if nextEnergy == 0.0
        delete!(PENDING_ATOM_DYNAMICS, atom)
        return nothing
    end
    dynamics = _pending_dynamics(atom)
    PENDING_ATOM_DYNAMICS[atom] = AtomDynamics(dynamics.velocityDirection, nextEnergy)
    return nothing
end

function GetNeighborVacancy(atom::Atom, simulator::Simulator)
    grid = simulator.grid    
    cell = GetCell(grid, atom.cellIndex)
    nearestVacancyDistance_squared = Inf
    nearestVacancyIndex = -1
    for neighborCellInfo in cell.neighborCellsInfo
        index, cross = neighborCellInfo.index, neighborCellInfo.cross
        neighborCell = GetCell(grid, index)
        for vacancy in neighborCell.vacancies
            dr2 = ComputeDistance_squared(atom.coordinate, vacancy.coordinate, cross, simulator.box)
            if dr2 < simulator.parameters.vacancyRecoverDistance_squared && dr2 < nearestVacancyDistance_squared
                nearestVacancyDistance_squared = dr2
                nearestVacancyIndex = vacancy.index
            end
        end
    end
    return nearestVacancyIndex
end


function Stop!(atom::Atom, simulator::Simulator)
    SetVelocityDirection!(atom, SVector{3,Float64}(0.0, 0.0, 0.0), simulator)
    SetEnergy!(atom, 0.0, simulator)
    Recover!(atom, simulator)
end


function Recover!(atom::Atom, simulator::Simulator)
    nearestVacancyIndex = GetNeighborVacancy(atom, simulator)
    if nearestVacancyIndex != -1
        SetOnLatticePoint!(atom, simulator.latticePoints[nearestVacancyIndex], simulator)
        deleteat!(simulator.vacancies, findfirst(v -> v.index == nearestVacancyIndex, simulator.vacancies))
        cell = GetCell(simulator.grid, atom.cellIndex)
        deleteat!(cell.vacancies, findfirst(v -> v.index == nearestVacancyIndex, cell.vacancies))
    end
end


function SetOnLatticePoint!(atom::Atom, latticePoint::LatticePoint, simulator::Simulator; isUpdateEnv::Bool = true)
    SetEnergy!(atom, 0.0)
    SetVelocityDirection!(atom, SVector{3,Float64}([0.0, 0.0, 0.0]))
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
    latticePoints = simulator.latticePoints
    for vacancy in simulator.vacancies
        latticePoint = latticePoints[vacancy.index]
        cell = GetCell(simulator.grid, latticePoint.cellIndex)
        empty!(cell.vacancies)
    end
    empty!(simulator.vacancies)

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
    dVectors = Vector{SVector{3,Float64}}()
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
    
    for i in eachindex(enviroment)
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
    if simulator.parameters.isAmorphous 
        rng = THREAD_RNG[Threads.threadid()]
        atom.coordinate .= GetCell(simulator.grid, atom.cellIndex).ranges[:,1] .+ [rand(rng) * simulator.grid.vectors[d, d] for d in 1:3]
    else
        ah = simulator.parameters.amorphousHeight
        if atom.coordinate[3] > ah
            rng = THREAD_RNG[Threads.threadid()]
            cell = GetCell(simulator.grid, atom.cellIndex)
            atom.coordinate[1] = cell.ranges[1,1] + rand(rng) * simulator.grid.vectors[1, 1]
            atom.coordinate[2] = cell.ranges[2,1] + rand(rng) * simulator.grid.vectors[2, 2]
            base = cell.ranges[3,1] > ah ? cell.ranges[3,1] : ah
            latticeTop = simulator.parameters.primaryVectors[3,3] * simulator.parameters.latticeRanges[3,2]     
            top = cell.ranges[3,2] < latticeTop ? cell.ranges[3,2] : latticeTop
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



function SetCoordinate!(atom::Atom, coordinate::AbstractVector{<:Real})
    if IS_DYNAMIC_LOAD
        atom.coordinate = SVector{3,Float64}(coordinate[1], coordinate[2], coordinate[3])
    else
        atom.coordinate .= coordinate
    end
end


function TemperatureToSigma(T::Float64, θ_D::Float64, m_rel::Float64; atol=1e-10, rtol=1e-8)
    if T == 0.0
        log_debug("Temperature is 0 K")
        return 0
    end
    ħ   = 1.054_571_817e-34      # J·s
    kB  = 1.380_649_000e-23      # J/K
    amu = 1.660_539_066_60e-27   # kg

    M = m_rel * amu
    y_max = θ_D / T      

    # Integration of x/(e^x-1) dx
    integrand(x) = x / (exp(x) - 1)
    I, _ = quadgk(integrand, 0.0, y_max; atol, rtol)

    σ2 = 3 * ħ^2 / (M * kB * θ_D) * (0.25 + (T/θ_D)^2 * I)
    σ  = sqrt(σ2) * 1e10         # m → Å

    return σ
end
