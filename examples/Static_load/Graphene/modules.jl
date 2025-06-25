function GetDTECustom(atom::Atom, simulator::Simulator)  # Modify. threshold displacement energy
    nlayers = 1
    b = 6.7
    z = atom.coordinate[2]
    isOnLatticeFlag = atom.latticePointIndex != -1 ? true : false
    if isOnLatticeFlag
        if z <= 2.25 * b || z > (1.25 + nlayers) * b 
            return 22.0
        elseif 2.25 * b < z <= 2.75 *b || (1+nlayers-0.25) * b < z <= (1+nlayers+0.25) * b 
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



function Irradiation(simulator::Simulator, energy::Float64)
    Restore!(simulator)
    simulator.nCascade += 1
    ionPosition = RandomInAnUnitGrapheneCell(1.42) + [18.855045, 22.981482313368623, 60]
    ion = Atom(3, ionPosition, parameters)
    SetVelocityDirection!(ion, [0.,0.,-1.])
    SetEnergy!(ion,energy)
    push!(simulator, ion)
    Cascade!(ion, simulator)    
    Vs, _ = DefectStatics(simulator)
    return length(Vs)
end