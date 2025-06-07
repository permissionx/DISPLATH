
function GetTargetsFromNeighbor_dynamicLoad_aborted(atom::Atom, gridCell::GridCell, simulator::Simulator)
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
            if neighborAtom.index == atom.index
                continue
            end
            Pertubation!(neighborAtom, simulator)
            if ComputeVDistance(atom, neighborAtom, neighborCellInfo.cross, box) > 0
                ComputeP!(atom, neighborAtom, neighborCellInfo.cross, box)
                if neighborAtom.pValue[atom.index] >= simulator.parameters.pMax
                    DeleteP!(neighborAtom, atom.index)
                    if atom.isNewlyLoaded == true
                        SetCoordinate!(atom, atom.latticeCoordinate)
                    end
                    continue
                end
                matchFlag = true
                for target in targets
                    if !SimultaneousCriteria(atom, neighborAtom, target, simulator)
                        matchFlag = false
                        DeleteP!(neighborAtom, atom.index)
                        if atom.isNewlyLoaded == true
                            SetCoordinate!(atom, atom.latticeCoordinate)
                        end
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
        Log("Infinitely fly atom in the $(simulator.nCascade)th irradiation:\n$(atom)\n")
    end
    return (targets, infiniteFlag)
end

function GetTargetsFromNeighbor_aborted(atom::Atom, gridCell::GridCell, simulator::Simulator)
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
        infiniteFlag = false
        for neighborAtom in neighborCell.atoms
            #if simulator.nCascade == 332 && neighborAtom.index == 296
            #    println(atom.index, " ", neighborAtom.index)
            #    println(neighborAtom)
            #    println(neighborCell.index)
            #    println(simulator.cellGrid.cells[11,21,11].atoms)
            #end
            if neighborAtom.index == atom.index
                continue
            end
            Pertubation!(neighborAtom, simulator)
            if ComputeVDistance(atom, neighborAtom, neighborCellInfo.cross, box) > 0 
                ComputeP!(atom, neighborAtom, neighborCellInfo.cross, box)
                if neighborAtom.pValue[atom.index] >= simulator.parameters.pMax
                    DeleteP!(neighborAtom, atom.index)
                    if neighborAtom.latticePointIndex != -1
                        SetCoordinate!(neighborAtom, simulator.latticePoints[neighborAtom.latticePointIndex].coordinate)
                    end
                    continue
                end
                matchFlag = true
                for target in targets
                    if !SimultaneousCriteria(atom, neighborAtom, target, simulator)
                        matchFlag = false
                        DeleteP!(neighborAtom, atom.index)
                        if neighborAtom.latticePointIndex != -1
                            SetCoordinate!(neighborAtom, simulator.latticePoints[neighborAtom.latticePointIndex].coordinate)
                        end
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
        Log("Infinitely fly atom in the $(simulator.nCascade)th irradiation:\n$(atom)\n")
    end
    return (targets, infiniteFlag)
end
