function GetCell(grid::Grid, cellIndex::Tuple{Int64, Int64, Int64})
    return _GetCellDense(grid, cellIndex)
end


function _GetCellDense(grid::Grid, cellIndex::Tuple{Int64, Int64, Int64})
    return grid.cells[cellIndex...]
end

function CreateCell(cellIndex::Tuple{Int64, Int64, Int64}, vectors::Matrix{Float64})
    x, y, z = cellIndex
    ranges = Matrix{Float64}(undef, 3, 2)
    ranges[1,1] = (x-1) * vectors[1,1]
    ranges[1,2] = x * vectors[1,1]
    ranges[2,1] = (y-1) * vectors[2,2]
    ranges[2,2] = y * vectors[2,2]
    ranges[3,1] = (z-1) * vectors[3,3]
    ranges[3,2] = z * vectors[3,3]  
    cell = Cell(cellIndex, Vector{Atom}(), Vector{LatticePoint}(), 
                            ranges, 
                            Array{NeighborCellInfo, 3}(undef, 3, 3, 3), false, 0.0)
    return cell
end

# belows are for dynamic load

function GetCell(grid::Grid, cellIndex::Tuple{Int64, Int64, Int64}, simulator::Simulator)
    return _GetCellDict!(grid, cellIndex, simualtor)
end

function _GetCellDict!(grid::Grid, cellIndex::Tuple{Int64, Int64, Int64}, simulator)
    dks = simulator.deprecatedCellKeys
    cells = grid.cells
    if haskey(cells, cellIndex)
        if cellIndex in dks
            delete!(dks, cellIndex)
        end
    else
        if isempty(dks)
            cells[cellIndex] = CreateCell(cellIndex, grid.vectors, simualtor)
        else
            dk = pop!(dks)
            cells[cellIndex] = pop!(cells, dk)
            UpdateCell!(cells[cellIndex])
        end
    end
    return cells[cellIndex]
end 

function InitCellStd!(simulator::Simulator, PN::Vector{Int64})
    # only for dyanmic load
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
                    atom.isLatticeAtom = true
                    push!(cellStd.atoms, atprimaryVectorsom)
                end
            end
        end
    end
    simulator.cellLatticeAtomNumber = indexInCell
end

function CreateCell(cellIndex::Tuple{Int64, Int64, Int64}, vectors::Matrix{Float64}, simulator::Simulator)
    parameters = simulator.parameters
    ranges = Matrix{Float64}(undef, 3, 2)
    for d in 1:3
        ranges[d,1] = (cellIndex[d] - 1) * vectors[1,1]
        ranges[d,2] = cellIndex[d] * vectors[1,1]
    end
    cell = Cell(cellIndex, Vector{Atom}(), Vector{LatticePoint}(), 
                            ranges, 
                            Array{NeighborCellInfo, 3}(undef, 3, 3, 3), false, 0.0)
    cell.isPushedNeighbor = false
    cell.hasNeighborObj = false
    for d in 1:3
        if ranges[d,1] < latticeRanges[d,1] * simulator.primaryVector[d,d] || ranges[d,1] > latticeRanges[d,2] * simulator.primaryVector[d,d]
            isEmpty = true
        end
    end
    for stdAtom in simulator.cellStd.atoms
        coords = [stdAtom.coordinate[d] + cell.ranges[d, 1] for d in 1:3]
        atom = Atom(stdAtom.type, coords, parameters)
        atom.indexInCell = stdAtom.indexInCell
        simulator.minLatticeAtomID -= 1
        atom.index = simulator.minLatticeAtomID
        atom.cellIndex = cellIndex
        atom.isLatticeAtom = true
        if isEmpty
            atom.isAlive = false
        else
            atom.isAlive = true
        end
        push(cell.latticeAtoms, atom)
    end
    return cell
end


function UpdateCell!(cell::Cell, cellIndex::Tuple{Int64, Int64, Int64}, vectors::Matrix{Float64}, simulator::Simulator)
    for d in 1:3
        cell.ranges[d, 1] = (cellIndex[d] - 1) * vectors[d,d]
        cell.ranges[d, 2] = cellIndex[d] * vectors[d,d]
    end
    cell.cellIndex = cellIndex
    cell.isPushedNeighbor = false  # update by create. 
    for d in 1:3
        if ranges[d,1] < latticeRanges[d,1] * simulator.primaryVector[d,d] || ranges[d,1] > latticeRanges[d,2] * simulator.primaryVector[d,d]
            isEmpty = true
        end
    end
    for atom, stdAtom in zip(cell.latticeAtoms, simulator.cellStd.atoms)
        [atom.coordinate[d] = stdAtom.coordinate[d] + cell.ranges[d,1] for d in 1:3]
        atom.cellIndex = cellIndex
        simulator.minLatticeAtomID -= 1
        atom.index = simulator.minLatticeAtomID
        if isEmpty
            atom.isAlive = false
        else
            atom.isAlive = true
        end
    end    
end



function GetNeighborCellsInfo!(cell, grid)
    if ! cell.isPushedNeighbor
        if ! cell.hasNeighborObj
            SetNeighborCellsInfo!(cell, grid)
        else
            UpdateNeighborCellsInfo!(cell, grid)
        end
    end
    return cell.neighborCellsInfo
end

function SetNeighborCellsInfo!(cell::Cell, grid::Grid)
    # Direct triple loop implementation - much faster than recursion
    for delta_x in [-1, 0, 1]
        for delta_y in [-1, 0, 1]
            for delta_z in [-1, 0, 1]
                neighborKeys = (Int8(delta_x), Int8(delta_y), Int8(delta_z))  
                neighborIndex = [0, 0, 0]  
                neighborCross = [Int8(0), Int8(0), Int8(0)]  
                # Calculate neighbor cell index and cross flags for each dimension
                for d in 1:3
                    delta = neighborKeys[d]
                    index = cell.index[d] + delta
                    cross = Int8(0)
                    if index < 1
                        index += grid.sizes[d]
                        cross = Int8(-1)
                    elseif index > grid.sizes[d]
                        index -= grid.sizes[d]
                        cross = Int8(1)
                    end
                    neighborIndex[d] = index
                    neighborCross[d] = cross
                end
                neighborIndex_tuple = (neighborIndex[1], neighborIndex[2], neighborIndex[3])
                neighborCross_tuple = (neighborCross[1], neighborCross[2], neighborCross[3])
                neighborCellInfo = NeighborCellInfo(neighborIndex_tuple, neighborCross_tuple)
                idx = (delta_x + 2, delta_y + 2, delta_z + 2)
                cell.neighborCellsInfo[idx...] = neighborCellInfo
            end
        end
    end
    cell.isPushedNeighbor = true 
    cell.hasNeighborObj = true
end



function UpdateNeighborCellsInfo!(cell::Cell, grid::Grid)
    # Direct triple loop implementation - much faster than recursion
    for delta_x in [-1, 0, 1]
        for delta_y in [-1, 0, 1]
            for delta_z in [-1, 0, 1]
                neighborCellInfo = cell.neighborCellsInfo[delta_x+2, delta_y+2, delta_z+2]
                neighborKeys = (Int8(delta_x), Int8(delta_y), Int8(delta_z))  
                for d in 1:3
                    delta = neighborKeys[d]
                    index = cell.index[d] + delta
                    cross = Int8(0)
                    if index < 1
                        index += grid.sizes[d]
                        cross = Int8(-1)
                    elseif index > grid.sizes[d]
                        index -= grid.sizes[d]
                        cross = Int8(1)
                    end
                    neighborCellInfo.index[d] = index
                    neighborCellInfo.cross[d] = cross
                end
            end
        end
    end
    cell.isPushedNeighbor = true 
end





function DeprecateCell!(cell::Cell, simulator::Simulator)
    if isempty(cell.atoms) && isempty(cell.vacancies)
        cell.isPushedNeighbor = false 
        delete!(simulator.deprecatedCellKeys, cell.index)
    end
end


function ChangeCell!(atom::Atom, nextCellIndex::Tuple{Int64, Int64, Int64}, simulator::Simulator)
    grid = simulator.grid
    if !IS_DYNAMIC_LOAD
        originalCell = GetCell(grid, atom.cellIndex, simulator)
        delete!(originalCell, atom, simulator)
        nextCell = GetCell(grid, nextCellIndex, simulator)
        push!(nextCell, atom, simulator)
    else
        originalCell = GetCell(grid, atom.cellIndex, simulator)
        delete!(originalCell, atom, simulator)
        nextCell = GetCell(grid, nextCellIndex, simulator)
        push!(nextCell, atom, simulator)
        nextCellIndexes = Set([neighborCellInfo.index for neighborCellInfo in GetNeighborCellsInfo!(nextCell, grid)])
        for cellInfo in GetNeighborCellsInfo!(originalCell, grid)
            neighborIndex = cellInfo.index
            if !neighborIndex in nextCellIndexes
                DeprecateCell!(GetCell(grid, neighborIndex, simulator), simulator)
            end
        end
            
    end
end