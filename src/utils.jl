function DefectStatics(simulator::Simulator)
    if !simulator.isStore
        error("simulator must be stored when counting defects")
    end
    latticePoints = simulator.latticePoints
    atoms = simulator.atoms
    vacancies = Vector{LatticePoint}()
    interstitials = Vector{Atom}()
    for idx in simulator.displacedAtoms
        if latticePoints[idx].atomIndex == -1
            push!(vacancies, latticePoints[idx])
        end
        if atoms[idx].isAlive && atoms[idx].latticePointIndex == -1
            push!(interstitials, atoms[idx])
        end
    end
    for idx in simulator.numberOfAtomsWhenStored:length(simulator.atoms)
        push!(interstitials, atoms[idx])
    end
    return vacancies, interstitials
end

function CountVacancies(simulator::Simulator)
    nVacancies = 0
    for latticePoint in simulator.latticePoints
        if latticePoint.atomIndex == -1
            nVacancies += 1
        end
    end
    return nVacancies
end

function ExtractVacancyLattices(simulator::Simulator)
    vacancyLattices = []
    for latticePoint in simulator.latticePoints
        if latticePoint.atomIndex == -1
            push!(vacancyLattices, latticePoint)
        end
    end
    return vacancyLattices
end


function RandomPointInCircle(radius::Float64=3.0)
    rng = THREAD_RNG[Threads.threadid()]
    θ = 2π * rand(rng)  
    r = sqrt(rand(rng)) * radius  
    x = r * cos(θ)
    y = r * sin(θ)
    return [x, y, 0.0]
end  

function RandomInSquare(a::Float64, b::Float64)
    rng = THREAD_RNG[Threads.threadid()]
    x = rand(rng) * a
    y = rand(rng) * b 
    return [x,y,0.0]
end

function RandomInAnUnitGrapheneCell(a::Float64)
    X = a * 3
    Y = sqrt(3) * a
    rng = THREAD_RNG[Threads.threadid()]
    x = rand(rng) * X
    y = rand(rng) * Y
    return [x, y, 0.0]
end  

function RandomVectorInUnitSphere(θrange::Float64)
    rng = THREAD_RNG[Threads.threadid()]
    θ = θrange * rand(rng)
    φ = 2π * rand(rng)
    x = sin(θ) * cos(φ)
    y = sin(θ) * sin(φ)
    z = cos(θ)
    return [x, y, -z]
end
