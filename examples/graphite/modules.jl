function GetDTECustom(atom::Atom, simulator::Simulator)  # Modify. threshold displacement energy
    b = 6.7
    z = atom.coordinate[2]
    isOnLatticeFlag = atom.latticePointIndex != -1 ? true : false
    if isOnLatticeFlag
        if z <= 2.25 * b || z > 4.25 * b 
            return 22.0
        elseif 2.25 * b < z <= 2.75 *b || 3.75 * b < z <= 4.25 * b 
            if atom.type == 1
                return 32.0
            else
                return 35.0
            end
        else
            if atom.type == 1
                return 42.0
            else
                return 45.0
            end
        end
    else
        return 2.0
    end
end

function RandomInAnUnitGrapheneCell(a::Float64)
    X = a * 3
    Y = sqrt(3) * a
    x = rand() * X
    y = rand() * Y
    return [x, y, 0.0]
end  

function CountVacancy(simulator::Simulator)
    latticePoints = simulator.latticePoints
    nV = 0
    for latticePoint in latticePoints
        if latticePoint.atomIndex == -1
            nV += 1
        end
    end
    return nV
end

function Irradiation(simulator::Simulator, energy::Float64)
    Restore!(simulator)
    simulator.nIrradiation += 1
    ionPosition = RandomInAnUnitGrapheneCell(1.42) + [18.855045, 22.981482313368623, 35]
    ion = Atom(3, ionPosition, parameters)
    SetVelocityDirection!(ion, [0.,0.,-1.])
    SetEnergy!(ion,energy)
    push!(simulator, ion)
    Cascade!(ion, simulator)    
    nV = CountVacancy(simulator)
    return nV
end