function GetCell(grid::Grid, cellIndex::Tuple{Int64, Int64, Int64})
    return _GetCellDense(grid, cellIndex)
end


function _GetCellDense(grid::Grid, cellIndex::Tuple{Int64, Int64, Int64})
    return grid.cells[cellIndex...]
end

function CreateCell(cellIndex::Tuple{Int64, Int64, Int64}, vectors::Matrix{Float64})
    ranges = Matrix{Float64}(undef, 3, 2)
    for d in 1:3
        ranges[d,1] = (cellIndex[d] - 1) * vectors[d,d]
        ranges[d,2] = cellIndex[d] * vectors[d,d]
    end
    cell = Cell(cellIndex, Vector{Atom}(), Vector{LatticePoint}(), 
                            ranges, 
                            Array{NeighborCellInfo, 3}(undef, 3, 3, 3), false, 0.0)
    return cell
end

# belows are for dynamic load

@inline function CellLower(cellIndex::Tuple{Int64, Int64, Int64}, d::Int64, grid::Grid)
    return (cellIndex[d] - 1) * grid.vectors[d,d]
end

@inline function CellUpper(cellIndex::Tuple{Int64, Int64, Int64}, d::Int64, grid::Grid)
    return cellIndex[d] * grid.vectors[d,d]
end

@inline CellLower(cell::Cell, d::Int64, grid::Grid) = CellLower(cell.index, d, grid)
@inline CellUpper(cell::Cell, d::Int64, grid::Grid) = CellUpper(cell.index, d, grid)

function CellRanges(cellIndex::Tuple{Int64, Int64, Int64}, grid::Grid)
    return @SMatrix [
        CellLower(cellIndex, 1, grid) CellUpper(cellIndex, 1, grid)
        CellLower(cellIndex, 2, grid) CellUpper(cellIndex, 2, grid)
        CellLower(cellIndex, 3, grid) CellUpper(cellIndex, 3, grid)
    ]
end

CellRanges(cell::Cell, grid::Grid) = CellRanges(cell.index, grid)

function IsEmptyDynamicCell(cellIndex::Tuple{Int64, Int64, Int64}, simulator::Simulator)
    latticeRanges = simulator.parameters.latticeRanges
    primaryVectors = simulator.parameters.primaryVectors
    grid = simulator.grid
    for d in 1:3
        lo = CellLower(cellIndex, d, grid)
        if lo < latticeRanges[d,1] * primaryVectors[d,d] || lo > latticeRanges[d,2] * primaryVectors[d,d]
            return true
        end
    end
    return false
end

@inline function CellKey(grid::Grid, cellIndex::Tuple{Int64, Int64, Int64})
    return ((cellIndex[1] - 1) * grid.sizes[2] + (cellIndex[2] - 1)) * grid.sizes[3] + cellIndex[3]
end

function GetCell(grid::Grid, cellIndex::Tuple{Int64, Int64, Int64}, simulator::Simulator)
    return _GetCellDict!(grid, cellIndex, simulator)
end

function _GetCellDict!(grid::Grid, cellIndex::Tuple{Int64, Int64, Int64}, simulator::Simulator)
    cells = grid.cells
    key = CellKey(grid, cellIndex)
    cell = get(cells, key, nothing)
    if cell === nothing
        if isempty(simulator.freeCells)
            cell = CreateCell(cellIndex, grid.vectors, simulator)
        else
            cell = pop!(simulator.freeCells)
            UpdateCell!(cell, cellIndex, grid.vectors, simulator)
        end
        cells[key] = cell
        push!(simulator.touchedCells, cell)
    end
    return cell
end

const FREE_CELL_POOL_MAX = 1 << 21

# Deferred reclamation: at cascade end, drop cells that are still empty from
# the grid dict and park them for reuse. Cells holding defects stay resident.
function SweepTouchedCells!(simulator::Simulator)
    grid = simulator.grid
    cells = grid.cells
    freeCells = simulator.freeCells
    for cell in simulator.touchedCells
        if isempty(cell.atoms) && isempty(cell.vacancies)
            key = CellKey(grid, cell.index)
            if get(cells, key, nothing) === cell
                delete!(cells, key)
                if length(freeCells) < FREE_CELL_POOL_MAX
                    push!(freeCells, cell)
                end
            end
        end
    end
    empty!(simulator.touchedCells)
end

# A cell from an earlier cascade can be emptied by the current one (its atom
# moved away or got deleted); record it so the sweep can reclaim it.
@inline function MarkCellIfEmpty!(cell::Cell, simulator::Simulator)
    if isempty(cell.atoms) && isempty(cell.vacancies)
        push!(simulator.touchedCells, cell)
    end
    return nothing
end

function InitCellStd!(simulator::Simulator, PN::Vector{Int64})
    # only for dyanmic load
    parameters = simulator.parameters
    basisTypes = parameters.basisTypes
    basis = parameters.basis
    primaryVectors = parameters.primaryVectors
    cellStd = simulator.cellStd
    indexInCell = 0
    cellStd = simulator.cellStd
    for X in 0:PN[1]-1
        for Y in 0:PN[2]-1
            for Z in 0:PN[3]-1
                for i in eachindex(basisTypes)
                    indexInCell += 1
                    x = (X + basis[i,1]) * primaryVectors[1,1]
                    y = (Y + basis[i,2]) * primaryVectors[2,2]
                    z = (Z + basis[i,3]) * primaryVectors[3,3]
                    atom = Atom(basisTypes[i], [x, y, z], parameters)
                    atom.index = 0 
                    push!(cellStd.atoms, atom)
                end
            end
        end
    end
    simulator.cellLatticeAtomNumber = indexInCell
end

function CreateCell(cellIndex::Tuple{Int64, Int64, Int64}, vectors::Matrix{Float64}, simulator::Simulator)
    return Cell(cellIndex, Atom[], Atom[])
end

function UpdateCell!(cell::Cell, cellIndex::Tuple{Int64, Int64, Int64}, vectors::Matrix{Float64}, simulator::Simulator)
    empty!(cell.atoms)
    empty!(cell.vacancies)
    cell.index = cellIndex
    return cell
end

function RefillLatticeAtoms!(cell::Cell, simulator::Simulator)
    return nothing
end

function GetNeighborCellsInfo!(cell, grid)
    if ! cell.isPushedNeighbor
        SetNeighborCellsInfo!(cell, grid)
    end
    return cell.neighborCellsInfo
end

function GetNeighborCellsInfo!(cell::Cell, grid::Grid, simulator, slot::Int=1)
    neighborCellsInfo = simulator.workBuffers.neighborCellsInfos[slot]
    UpdateNeighborCellsInfo!(neighborCellsInfo, cell.index, grid)
    return neighborCellsInfo
end

@inline function _neighbor_index_cross(index::Int64, delta::Int64, size::Int64)
    neighborIndex = index + delta
    cross = Int8(0)
    if neighborIndex < 1
        neighborIndex += size
        cross = Int8(-1)
    elseif neighborIndex > size
        neighborIndex -= size
        cross = Int8(1)
    end
    return neighborIndex, cross
end

function SetNeighborCellsInfo!(cell::Cell, grid::Grid)
    # Direct triple loop implementation - much faster than recursion
    for delta_x in -1:1
        for delta_y in -1:1
            for delta_z in -1:1
                ix, cx = _neighbor_index_cross(cell.index[1], delta_x, grid.sizes[1])
                iy, cy = _neighbor_index_cross(cell.index[2], delta_y, grid.sizes[2])
                iz, cz = _neighbor_index_cross(cell.index[3], delta_z, grid.sizes[3])
                neighborCellInfo = NeighborCellInfo((ix, iy, iz), (cx, cy, cz))
                idx = (delta_x + 2, delta_y + 2, delta_z + 2)
                cell.neighborCellsInfo[idx...] = neighborCellInfo
            end
        end
    end
    cell.isPushedNeighbor = true 
end

function UpdateNeighborCellsInfo!(
    neighborCellsInfo::Array{NeighborCellInfo, 3},
    cellIndex::Tuple{Int64, Int64, Int64},
    grid::Grid,
)
    for delta_x in -1:1
        for delta_y in -1:1
            for delta_z in -1:1
                neighborCellInfo = neighborCellsInfo[delta_x+2, delta_y+2, delta_z+2]
                ix, cx = _neighbor_index_cross(cellIndex[1], delta_x, grid.sizes[1])
                iy, cy = _neighbor_index_cross(cellIndex[2], delta_y, grid.sizes[2])
                iz, cz = _neighbor_index_cross(cellIndex[3], delta_z, grid.sizes[3])
                neighborCellInfo.index = (ix, iy, iz)
                neighborCellInfo.cross = (cx, cy, cz)
            end
        end
    end
    return neighborCellsInfo
end



function UpdateNeighborCellsInfo!(cell::Cell, grid::Grid)
    # Direct triple loop implementation - much faster than recursion
    for delta_x in -1:1
        for delta_y in -1:1
            for delta_z in -1:1
                neighborCellInfo = cell.neighborCellsInfo[delta_x+2, delta_y+2, delta_z+2]
                ix, cx = _neighbor_index_cross(cell.index[1], delta_x, grid.sizes[1])
                iy, cy = _neighbor_index_cross(cell.index[2], delta_y, grid.sizes[2])
                iz, cz = _neighbor_index_cross(cell.index[3], delta_z, grid.sizes[3])
                neighborCellInfo.index = (ix, iy, iz)
                neighborCellInfo.cross = (cx, cy, cz)
            end
        end
    end
    cell.isPushedNeighbor = true 
end




function ChangeCell!(atom::Atom, nextCellIndex::Tuple{Int64, Int64, Int64}, simulator::Simulator)
    grid = simulator.grid
    if !IS_DYNAMIC_LOAD
        originalCell = GetCell(grid, atom.cellIndex)
        delete!(originalCell, atom, simulator)
        nextCell = GetCell(grid, nextCellIndex)
        push!(nextCell, atom, simulator)
    else
        originalCell = GetCell(grid, atom.cellIndex, simulator)
        delete!(originalCell, atom, simulator)
        MarkCellIfEmpty!(originalCell, simulator)
        nextCell = GetCell(grid, nextCellIndex, simulator)
        push!(nextCell, atom, simulator)
    end
end
