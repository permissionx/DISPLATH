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

function AddKMCInfo(atom::Atom, simulator::Simulator)
    if atom.latticePointIndex == -1
        error("Atom is not on lattice")
    end
    latticePoint = simulator.latticePoints[atom.latticePointIndex]
    envIndex = GetEnviromentIndex(latticePoint, simulator)
    finalLatticePointEnvIndexs, finalLatticePointIndexs = GetFinalLatticePointInfo(latticePoint, simulator)
    energyBarriers = [simulator.energyBarrierDict[envIndex, i] for i in finalLatticePointEnvIndexs]
    
    atom.probabilities = [ComputeProbability(e, simulator) for e in energyBarriers]
    atom.probability = sum(atom.probabilities)
    atom.finalLatticePointIndexs = finalLatticePointIndexs

    simulator.probability += atom.probability
    push!(simulator.probabilities, atom.probability)
    push!(simulator.mobileAtoms, atom)
end 

function ComputeProbability(energyBarrier::Float64, simulator::Simulator)
    return exp(-energyBarrier / simulator.kmc.temperature_kb)
end

function GetRandomMobileAtom(simulator::Simulator)
    randomNumber = rand() * simulator.probability
    cumulativeProbability = 0.0
    for i in 1:length(simulator.probabilities)
        cumulativeProbability += simulator.probabilities[i]
        if randomNumber <= cumulativeProbability
            return simulator.mobileAtoms[i], i
        end
    end
end

function GetRandomFinalLatticePointIndex(atom::Atom)
    randomNumber = rand() * atom.probability
    cumulativeProbability = 0.0
    for i in 1:length(atom.probabilities)
        cumulativeProbability += atom.probabilities[i]
        if randomNumber <= cumulativeProbability
            return atom.finalLatticePointIndexs[i]
        end
    end
end

function RandomMobileEvent(simulator::Simulator)
    atom, mobileIndex = GetRandomMobileAtom(simulator)
    latticePoint = simulator.latticePoints[GetRandomFinalLatticePointIndex(atom)]
    SetOnLatticePoint!(atom, latticePoint, simulator)
    dt = GetTimeStep(simulator)
    simulator.time += dt    
    DeleteKMCInfo!(simulator, mobileIndex)

    envIndexs = GetEnviromentIndex(latticePoint, simulator)
    for index in envIndexs
        latticePoint = simulator.latticePoints[index]
        atom = simulator.atoms[latticePoint.atomIndex]
        AddKMCInfo(atom, simulator)
    end
end

function DeleteKMCInfo!(simulator::Simulator, index::Int64)
    simulator.probability -= simulator.probabilities[index]
    deleteat!(simulator.probabilities, index)
    deleteat!(simulator.mobileAtoms, index)
end 

function UpdateKMCEnvironment(latticePoint::LatticePoint, simulator::Simulator)
     for envLatticePoint in latticePoint
        if envLatticePoint.atomIndex != -1
            atom = simulator.atoms[envLatticePoint.atomIndex]
            if atom.eventIndex != -1
                UpdateKMCInfo(atom, simulator)
            else
                AddKMCInfo(atom, simulator)
            end
        end
    end
end

    

function UpdateAtomKMCInfo(atom::Atom, simulator::Simulator)
    latticePoint = simulator.latticePoints[atom.latticePointIndex]
    envIndex = GetEnviromentIndex(latticePoint, simulator)
    if envIndex == simulator.perfectEnvIndex
        DeleteKMCInfo(simulator, atom.envIndex)
    else
        simulator.probabilities -= atom.probability
        finalLatticePointEnvIndexs, finalLatticePointIndexs = GetFinalLatticePointInfo(latticePoint, simulator)
        energyBarriers = [simulator.energyBarrierDict[envIndex, i] for i in finalLatticePointEnvIndexs]
        atom.probabilities = [ComputeProbability(e, simulator) for e in energyBarriers]
        atom.probability = sum(atom.probabilities)
        atom.finalLatticePointIndexs = finalLatticePointIndexs


x

# types: simualtor.kmc, simulator.kmc.probabilities,simulator.kmc.finalLatticePoints, simulator.kmc.mobileAtomindexs
# functions: GetTimeStep  

