function UpdateEvents!(latticePointIndexs::Set{Int64}, simulator::Simulator)
    if simulator.parameters.isKMC
        for latticePointIndex in latticePointIndexs
            latticePoint = simulator.latticePoints[latticePointIndex]
            if latticePoint.atomIndex != -1
                envIndex = GetEnvironmentIndex(latticePoint, simulator)
                atom = simulator.atoms[latticePoint.atomIndex]
                if envIndex == simulator.parameters.perfectEnvIndex
                    DeleteAtomEvents!(simulator, atom)
                else
                    if atom.eventIndex != -1
                        simulator.frequency -= atom.frequency
                        SetAtomEvents!(atom, latticePoint, simulator)
                        simulator.frequency += atom.frequency
                        simulator.frequencies[atom.eventIndex] = atom.frequency
                    else
                        SetAtomEvents!(atom, latticePoint, simulator)
                        AppendAtomEvents!(simulator, atom)
                    end
                end
            end
        end
    end
end


function SetAtomEvents!(atom::Atom, latticePoint::LatticePoint, simulator::Simulator)
    finalLatticePointEnvIndexs, finalLatticePointIndexs = GetFinalLatticePointInfo(latticePoint, simulator)
    energyBarriers = [simulator.energyBarrierDict[atom.type][[envIndex, i]] for i in finalLatticePointEnvIndexs]
    atom.frequencies = [ComputeKMCProbability(type, e, simulator.temperature) for e in energyBarriers]
    atom.frequency = sum(atom.frequencies)
    atom.finalLatticePointIndexs = finalLatticePointIndexs
end


function AppendAtomEvents!(simulator::Simulator, atom::Atom)
    push!(simulator.frequencies, atom.frequency)
    push!(simulator.mobileAtoms, atom)
    simulator.maxEventIndex += 1
    atom.eventIndex = simulator.maxEventIndex 
    simulator.frequency += atom.frequency
end 


function GetFinalLatticePointInfo(latticePoint::LatticePoint, simulator::Simulator)
    environment = latticePoint.environment
    latticePoints = simulator.latticePoints
    finalLatticePointEnvIndexs = Vector{Int64}()
    finalLatticePointIndexs = Vector{Int64}()
    for i in 1:length(environment)
        neighborLatticePoint = latticePoints[environment[i]]
        if neighborLatticePoint.atomIndex == -1 &&  latticePoint.tpye == atom.type
            push!(finalLatticePointEnvIndexs, 2^(i-1))    
            push!(finalLatticePointIndexs, environment[i])
        end
    end
    return finalLatticePointEnvIndexs, finalLatticePointIndexs
end

function SetTemperature!(simulator::Simulator, temperature::Float64)
    temperature_kb = temperature * 8.61733362E-5
    simulator.parameters.temperature_kb = temperature_kb
end

function ComputeProbability(type::Int64, energyBarrier::Float64, simulator::Simulator)
    return simulator.parameters.nu_0_dict[type] * exp(-energyBarrier / simulator.parameters.temperature_kb)
end

function GetRandomMobileAtomIndex(simulator::Simulator)
    randomNumber = rand() * (simulator.frequency + simulator.parameters.irrdiationFrequency)
    cumulativeProbability = 0.0
    for i in 1:length(simulator.frequencies)
        cumulativeProbability += simulator.frequencies[i]
        if cumulativeProbability >= randomNumber
            return i
        end
    end
    return -1
end

function GetRandomFinalLatticePointIndex(atom::Atom)
    randomNumber = rand() * atom.frequency
    cumulativeProbability = 0.0
    for i in 1:length(atom.frequencies)
        cumulativeProbability += atom.frequencies[i]
        if randomNumber <= cumulativeProbability
            return atom.finalLatticePointIndexs[i]
        end
    end
end

function ProcessAnEvent(simulator::Simulator)
    idx = GetRandomMobileAtomIndex(simulator)
    if idx > 0
        atom = simulator.mobileAtoms[idx]
        latticePoint = simulator.latticePoints[GetRandomFinalLatticePointIndex(atom)]
        Migrate!(atom, latticePoint, simulator)
    else
        Irrdiate!(simulator)
    end
    ElapsTime(simulator)
end


function DeleteAtomEvents!(simulator::Simulator, atom::Atom)
    index = atom.eventIndex
    simulator.frequency -= simulator.frequencies[index]
    deleteat!(simulator.frequencies, index)
    deleteat!(simulator.mobileAtoms, index)
    atom.eventIndex = -1
end 


function ElapsTime(simulator::Simulator)
    time = -1/simulator.frequency * log(rand())
    simulator.time += time
end


function Migrate!(atom::Atom, latticePoint::LatticePoint, simulator::Simulator)
    oldLatticePoint = simulator.latticePoints[atom.latticePointIndex]
    LeaveLatticePoint!(atom, simulator; isUpdateEnv = false)
    SetOnLatticePoint!(atom, latticePoint, simulator; isUpdateEnv = false)
    latticePointIndexs = Set([oldLatticePoint.environment; latticePoint.environment; latticePoint.index])
    UpdateEvents!(latticePointIndexs, simulator)
end



