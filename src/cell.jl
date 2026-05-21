function GetCell(grid::Grid, cellIndex::Tuple{Int64, Int64, Int64})
    if !IS_DYNAMIC_LOAD
        return _GetCellDense(grid, cellIndex)
    else
        return _GetCellDict!(grid, cellIndex)
    end
end

function _GetCellDense(grid::Grid, cellIndex::Tuple{Int64, Int64, Int64})
    return grid.cells[cellIndex...]
end

# belows are only for dynamic load 
function _GetCellDict!(grid::Grid, cellIndex::Tuple{Int64, Int64, Int64})
    dks = simulator.deprecatedCellKeys
    cells = grid.cells
    if haskey(cells, cellIndex)
        if cellIndex in dks
            delete!(dks, cellIndex)
        end
    else
        if isempty(dks)
            cells[cellIndex] = CreateCell(cellIndex, grid.vectors)
        else
            dk = pop!(dks)
            cells[cellIndex] = pop!(cells, dk)
            UpdateCell!(cells[cellIndex])
        end
    end
    return cells[cellIndex]
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
    cell.hasNeighbor = false
    for d in 1:3
        if ranges[d,1] < latticeRanges[d,1] * simulator.primaryVector[d,d] || ranges[d,1] > latticeRanges[d,2] * simulator.primaryVector[d,d]
            isEmpty = true
        end
    end
    for stdAtom in simulator.cellStd.atoms
        coords = [stdAtom.coordinate[d] + cell.ranges[d, 1] for d in 1:3]
        atom = Atom(stdAtom.type, coords, parameters)
        atom.indexInCell = stdAtom.indexInCell
        atom.index = 0
        atom.cellIndex = cellIndex
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
        if isEmpty
            atom.isAlive = false
        else
            atom.isAlive = true
        end
    end    
end


function SetCellNeighborInfo!(cell::Cell, grid::Grid)
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
    cell.hasNeighbor = true
end

function UpdateCellNeighborInfo!(cell::Cell, grid::Grid)
    # Direct triple loop implementation - much faster than recursion
    for delta_x in [-1, 0, 1]
        for delta_y in [-1, 0, 1]
            for delta_z in [-1, 0, 1]
                neighborCellsInfo = cell.neighborCellsInfo[delta_x+2, delta_y+2, delta_z+2]
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


function CreateGrid(box::Box, inputVectors::Matrix{Float64})
    if !box.isOrthogonal
        error("The box is not orthogonal, please use the orthogonal box.")
    end
    sizes = Vector{Int64}(undef, 3)
    vectors = zeros(Float64, 3, 3)
    for d in 1:3
        sizes[d] = Int64(floor(box.vectors[d,d] / inputVectors[d,d]))
        if sizes[d] < 3
            error("The box size in dimension $d is too small,  use a larger box! (At least 3 cells in each dimension)")
            exit()
        end
        vectors[d,d] = box.vectors[d,d] / sizes[d]
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
            SetCellNeighborInfo!(cell, grid)
        end
    else
        cells = Dict{Tuple{Int64, Int64, Int64}, Cell}()    
        cellVolume = vectors[1,1] * vectors[2,2] * vectors[3,3]
        grid = Grid(cells, vectors, sizes, cellVolume) 
    end
    log_success("Cell grid created")
    log_separator()
    return grid
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
                    push!(cellStd.atoms, atprimaryVectorsom)
                end
            end
        end
    end
    simulator.cellLatticeAtomNumber = indexInCell
end

function DeprecatCell!(cell::Cell, simulator::Simulator)
    if isEmpty(cell.atoms) && isEmpty(cell.vacancies)
        cell.isPushedNeighbor = false 
        delete!(simulator.deprecatedCellKeys, cell.index)
    end
end

