using PyCall
@pyimport ase
@pyimport dscribe.descriptors as descriptors

function GetDTE(atom::Atom, simulator::Simulator)
    if simulator.parameters.DTEMode == 1  # direct 
        return atom.dte
    elseif simulator.parameters.DTEMode == 2  # all environment
        return GetDTEByEnvironment(atom, simulator)
    elseif simulator.parameters.DTEMode == 3   # soap
        return GetDTEBySoap(atom, simulator)
    end
end

function GetBDE(atom::Atom, simulator::Simulator)  # BDE: binding energy
    if simulator.parameters.DTEMode == 1  # direct 
        return atom.bde
    elseif simulator.parameters.DTEMode == 2  # all environment
        return GetBDEByEnvironment(atom, simulator)
    elseif simulator.parameters.DTEMode == 3   # soap
        return GetBDEBySoap(atom, simulator)
    end
end

  
function GetDTE(latticePoint::LatticePoint, simulator::Simulator)
    index = GetEnvironmentIndex(latticePoint, simulator) 
    return simulator.DTEData[latticePoint.type][index]
end

function GetDTEByEnvironment(atom::Atom, simulator::Simulator)\
    if atom.latticePointIndex != -1 
        return GetDTE(simulator.latticePoints[atom.latticePointIndex], simulator)
    else
        return atom.dte/2
    end
end
         
function GetBDEByEnvironment(atom::Atom, simulator::Simulator)
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

