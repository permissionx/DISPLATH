function DefectStatics(simulator::Simulator)
    if !IS_DYNAMIC_LOAD
        if !simulator.isStore
            error("simulator must be stored when counting defects")
        end
        latticePoints = simulator.latticePoints
        atoms = simulator.atoms
        vacancies = Vector{LatticePoint}()
        interstitials = Vector{Atom}()
        for idx in simulator.displacedAtoms
            if latticePoints[idx].atomIndex == -1 || atoms[latticePoints[idx].atomIndex].type == length(keys(simulator.parameters.typeDict))
                push!(vacancies, latticePoints[idx])
            end
            if atoms[idx].isAlive && atoms[idx].latticePointIndex == -1
                push!(interstitials, atoms[idx])
            end
        end
        for idx in simulator.numberOfAtomsWhenStored:simulator.maxAtomID
            if atoms[idx].isAlive && atoms[idx].latticePointIndex == -1
                push!(interstitials, atoms[idx])
            end
        end
    else
        interstitials, vacancies = (simulator.atoms, simulator.vacancies)
    end
    return interstitials, vacancies
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


function rotation_matrix_from_vectors(vec1::AbstractVector, vec2::AbstractVector)
    a = normalize(vec1)
    b = normalize(vec2)
    v = cross(a, b)
    c = dot(a, b)

    if c ≈ 1.0
        return I(3) 
    end


    if c ≈ -1.0

        other = abs(dot(a, [0, 0, 1])) < 0.9 ? [0, 0, 1] : [1, 0, 0]
        axis = normalize(cross(a, other))
        return AngleAxis(π, axis...)|>RotMatrix
    end
    

    s = norm(v)
    vx = [0 -v[3] v[2]; v[3] 0 -v[1]; -v[2] v[1] 0]
    R = I(3) + vx + vx^2 * (1 / (1 + c))
    
    return R
end


function RandomlyDeviatedVector(incident_direction::AbstractVector, θrange::Float64)
    rng = THREAD_RNG[Threads.threadid()]
    θ = θrange * rand(rng)
    φ = 2π * rand(rng)
    
    x = sin(θ) * cos(φ)
    y = sin(θ) * sin(φ)
    z = cos(θ)
    

    local_vec = [x, y, -z]


    standard_axis = [0.0, 0.0, -1.0]
    R = rotation_matrix_from_vectors(standard_axis, incident_direction)

    final_vec = R * local_vec
    
    return final_vec
end