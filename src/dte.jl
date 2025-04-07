using PyCall
@pyimport ase
@pyimport dscribe.descriptors as descriptors

function GetDTE(atom::Atom, simulator::Simulator)
    if simulator.DTEMode == 1  # direct 
        return atom.dte
    elseif simulator.DTEMode == 2  # all enviroment
        return GetDTEByEnviroment(atom, simulator)
    elseif simulator.DTEMode == 3   # soap
        return GetDTEBySoap(atom, simulator)
    end
end

function GetBDE(atom::Atom, simulator::Simulator)  # BDE: binding energy
    if simulator.dteMode == 1  # direct 
        return atom.bde
    elseif simulator.DTEMode == 2  # all enviroment
        return GetBDEByEnviroment(atom, simulator)
    elseif simulator.DTEMode == 3   # soap
        return GetBDEBySoap(atom, simulator)
    end
end

function GetEnviromentLatticePoints(latticePoint::LatticePoint, simulator::Simulator)
    cellIndex = latticePoint.cellIndex
    theCell = simulator.cellGrid.cells[cellIndex[1], cellIndex[2], cellIndex[3]]
    cells = simulator.cellGrid.cells
    cut_squared = simulator.enviromentCut^2
    box = simulator.box
    enviromentLatticePointsIndex = Vector{Int64}()
    dVectors = Vector{Vector{Float64}}()
    allLatticePoints = simulator.latticePoints
    for (_, neighborCellInfo) in theCell.neighborCellsInfo
        index, cross = neighborCellInfo.index, neighborCellInfo.cross
        cell = cells[index[1], index[2], index[3]]
        latticePointsIndex = cell.latticePoints
        for neighborLatticePointIndex in latticePointsIndex
            neighborLatticePoint = allLatticePoints[neighborLatticePointIndex]
            if ComputeDistance_squared(latticePoint.coordinate, neighborLatticePoint.coordinate, cross, box) <= cut_squared && neighborLatticePointIndex != latticePoint.index
                push!(dVectors, VectorDifference(latticePoint.coordinate, neighborLatticePoint.coordinate, cross, box))
                push!(enviromentLatticePointsIndex, neighborLatticePointIndex)
            end
        end
    end
    # sort indexes by x then y then z of dVectors
    sorted_indices = sortperm(dVectors, by = v -> (v[1], v[2], v[3]))
    enviromentLatticePointsIndex = enviromentLatticePointsIndex[sorted_indices]
    
    return enviromentLatticePointsIndex
end

function InitLatticePointEnvronment(simulator::Simulator)
    for latticePoint in simulator.latticePoints
        latticePoint.enviroment = GetEnviromentLatticePoints(latticePoint, simulator)
    end
    # simulator.enviromentLength should be get from the DTEDict.
end


function GetEnviromentIndex(latticePoint::LatticePoint, simulator::Simulator)
    enviroment = latticePoint.enviroment
    latticePoints = simulator.latticePoints
    index = 0
    
    for i in 1:length(enviroment)
        if latticePoints[enviroment[i]].atomIndex != -1
            index += 2^(i-1)
        end
    end
    return index 
end 
  
function GetDTE(latticePoint::LatticePoint, simulator::Simulator)
    index = GetEnviromentIndex(latticePoint, simulator)
    return simulator.DTEData[latticePoint.type][index]
end

function GetDTEByEnviroment(atom::Atom, simulator::Simulator)\
    try        
        if atom.latticePointIndex != -1 
            return GetDTE(simulator.latticePoints[atom.latticePointIndex], simulator)
        else
            return atom.dte/2
        end
    catch
        return atom.dte
    end
end
       
function GetBDEByEnviroment(atom::Atom, simulator::Simulator)
    return GetDTE(atom, simulator)
end



function InitSoap(parameters::Parameters)
    soapParameters = parameters.soapParameters 
    soap = descriptors.SOAP(
        species = [type.name for type in values(parameters.typeDict)],
        periodic = true,
        r_cut = soapParameters[1],
        n_max = Int(soapParameters[2]),
        l_max = Int(soapParameters[3])
    )
    return soap
end

function GetNeighborArray(atom::Atom, simulator::Simulator)
    coordinates = atom.coordinate'  
    elementNames = Vector{String}([simulator.parameters.typeDict[atom.type].name])
    atoms = simulator.atoms
    typeDict = simulator.parameters.typeDict
    cellGrid = simulator.cellGrid
    theCell = cellGrid.cells[atom.cellIndex[1], atom.cellIndex[2], atom.cellIndex[3]]
    for (_, neighborInfo) in theCell.neighborCellsInfo
        index = neighborInfo.index
        cell = cellGrid.cells[index[1], index[2], index[3]]
        for atomIndex in cell.atoms
            if atomIndex != atom.index
                coordinate = atoms[atomIndex].coordinate'
                coordinates = vcat(coordinates, coordinate)
                elementNames = push!(elementNames, typeDict[atoms[atomIndex].type].name)
            end
        end
    end
    return coordinates, elementNames
end

       
function CreateSoap(atom::Atom, simulator::Simulator)
    cell = [simulator.box.vectors[i,i] for i in 1:3]
    coordinates, elementNames = GetNeighborArray(atom, simulator)
    pbc = true
    atoms = ase.Atoms(elementNames, coordinates, cell = cell, pbc = pbc)
    soap = simulator.soap

    descriptor = soap.create(atoms, centers=[0])
    return soap, descriptor
end

