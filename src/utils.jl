function DefectStatics(simulator::Simulator)
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
    θ = 2π * rand()  
    r = sqrt(rand()) * radius  
    x = r * cos(θ)
    y = r * sin(θ)
    return [x, y, 0.0]
end  

function RandomInSquare(a::Float64, b::Float64)
    x = rand() * a
    y = rand() * b 
    return [x,y,0.0]
end

function RandomInAnUnitGrapheneCell(a::Float64)
    X = a * 3
    Y = sqrt(3) * a
    x = rand() * X
    y = rand() * Y
    return [x, y, 0.0]
end  

function RandomVectorInUnitSphere(θrange::Float64)
    θ = θrange * rand()
    φ = 2π * rand()
    x = sin(θ) * cos(φ)
    y = sin(θ) * sin(φ)
    z = cos(θ)
    return [x, y, -z]
end
