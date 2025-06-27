#using PyCall
#@pyimport ase
#@pyimport dscribe.descriptors as descriptors
#
function GetDTE(atom::Atom, simulator::Simulator)
    if simulator.parameters.DTEMode == 1  # direct 
        return atom.dte
    elseif simulator.parameters.DTEMode == 2  # all environment
        return GetDTEByEnvironment(atom, simulator)
    #elseif simulator.parameters.DTEMode == 3   # soap
    #    return GetDTEBySoap(atom, simulator)
    elseif simulator.parameters.DTEMode == 4
        return GetDTECustom(atom, simulator)
    end
end

function GetBDE(atom::Atom, simulator::Simulator)  # BDE: binding energy
    if simulator.parameters.DTEMode == 1  # direct 
        return atom.bde
    elseif simulator.parameters.DTEMode == 2  # all environment
        return GetBDEByEnvironment(atom, simulator)
    #elseif simulator.parameters.DTEMode == 3   # soap
    #    return GetBDEBySoap(atom, simulator)
    elseif simulator.parameters.DTEMode == 4
        return GetBDECustom(atom, simulator)
    end
end

  
function GetDTE(latticePoint::LatticePoint, simulator::Simulator)
    index = GetEnvironmentIndex(latticePoint, simulator) 
    return simulator.DTEData[latticePoint.type][index]
end

function GetDTEByEnvironment(atom::Atom, simulator::Simulator)\
    if atom.latticePointIndex != -1 && simulator.latticePoints[atom.latticePointIndex].type == atom.type
        return GetDTE(simulator.latticePoints[atom.latticePointIndex], simulator)
    else
        return 0.1
    end
end
         
function GetBDEByEnvironment(atom::Atom, simulator::Simulator)
    return GetDTE(atom, simulator)
end

function GetBDECustom(atom::Atom, simulator::Simulator)
    return GetDTECustom(atom, simulator)
end


#function InitSoap(parameters::Parameters)
#    soapParameters = parameters.soapParameters 
#    soap = descriptors.SOAP(
#        species = [type.name for type in values(parameters.typeDict)],
#        periodic = true,
#        r_cut = soapParameters[1],
#        n_max = Int(soapParameters[2]),
#        l_max = Int(soapParameters[3])
#    )
#    return soap
#end
#test
function GetNeighborArray(atom::Atom, simulator::Simulator)
    coordinates = atom.coordinate'  
    elementNames = Vector{String}([simulator.parameters.typeDict[atom.type].name])
    typeDict = simulator.parameters.typeDict
    cellGrid = simulator.cellGrid
    theCell = cellGrid.cells[atom.cellIndex...]
    for neighborInfo in theCell.neighborCellsInfo
        index = neighborInfo.index
        cell = cellGrid.cells[index...]
        for cellAtom in cell.atoms
            if cellAtom.index != atom.index
                coordinate = cellAtom.coordinate'
                coordinates = vcat(coordinates, coordinate)
                elementNames = push!(elementNames, typeDict[cellAtom.type].name)
            end
        end
    end
    return coordinates, elementNames
end

       
#function CreateSoap(atom::Atom, simulator::Simulator)
#    cell = [simulator.box.vectors[i,i] for i in 1:3]
#    coordinates, elementNames = GetNeighborArray(atom, simulator)
#    pbc = true
#    atoms = ase.Atoms(elementNames, coordinates, cell = cell, pbc = pbc)
#    soap = simulator.soap
#
#    descriptor = soap.create(atoms, centers=[0])
#    return soap, descriptor
#end

