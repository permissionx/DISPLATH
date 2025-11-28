function Dump(simulator::Simulator, fileName::String, step::Int64, type::String="a", isDebug::Bool=false)
    if !simulator.parameters.isOrthogonal        
        error("The box is not orthogonal, please use the orthogonal box.")
    end
    open(fileName, type) do file
        write(file, "ITEM: TIMESTEP\n")
        write(file, string(step), "\n")
        write(file, "ITEM: NUMBER OF ATOMS\n")
        write(file, string(simulator.numberOfAtoms), "\n")
        write(file, "ITEM: BOX BOUNDS ")
        for d in 1:3
            if simulator.parameters.periodic[d] 
                write(file, "pp ")
            else
                write(file, "ff ")
            end
        end
        write(file, "\n")
        for d in 1:3
            write(file, "0 $(simulator.box.vectors[d,d])\n")
        end
        if isDebug
            write(file, "ITEM: ATOMS id type x y z vx vy vz energy cx cy cz dte\n")
        else
            write(file, "ITEM: ATOMS id type x y z e\n")
        end
        for atom in simulator.atoms
            if atom.isAlive
                if isDebug
                    write(file, "$(atom.index) $(atom.type) \
                    $(atom.coordinate[1]) $(atom.coordinate[2]) $(atom.coordinate[3]) \
                    $(atom.velocityDirection[1]*sqrt(2*atom.mass*atom.energy)) $(atom.velocityDirection[2]*sqrt(2*atom.mass*atom.energy)) $(atom.velocityDirection[3]*sqrt(2*atom.mass*atom.energy)) \
                    $(atom.energy) \
                    $(atom.cellIndex[1]) $(atom.cellIndex[2]) $(atom.cellIndex[3]) \
                    $(GetDTE(atom, simulator))\n")
                else
                    write(file, "$(atom.index) $(atom.type) \
                    $(atom.coordinate[1]) $(atom.coordinate[2]) $(atom.coordinate[3]) $(atom.energy)\n")
                end
            end
        end 
    end
end

function ReadDate(fileName::String, replicate::Vector{Int64})
    if length(replicate) != 3
        error("Replicate vector must contain 3 entries for x/y/z directions.")
    end
    if any(r -> r < 1, replicate)
        error("Replicate factors must be >= 1.")
    end
    xlo = NaN; xhi = NaN
    ylo = NaN; yhi = NaN
    zlo = NaN; zhi = NaN
    vectors = zeros(Float64, 3, 3)
    hasVector = falses(3)
    origin = zeros(Float64, 3)
    hasOrigin = false
    types = Int64[]
    xs = Float64[]
    ys = Float64[]
    zs = Float64[]
    open(fileName, "r") do f
        lines = readlines(f)
        i = 1
        while i <= length(lines)
            line = strip(lines[i])
            if isempty(line) || startswith(line, "#")
                i += 1
                continue
            end
            words = split(line)
            if length(words) >= 4 && words[4] == "xhi"
                xlo = parse(Float64, words[1])
                xhi = parse(Float64, words[2])
            elseif length(words) >= 4 && words[4] == "yhi"
                ylo = parse(Float64, words[1])
                yhi = parse(Float64, words[2])
            elseif length(words) >= 4 && words[4] == "zhi"
                zlo = parse(Float64, words[1])
                zhi = parse(Float64, words[2])
            elseif length(words) >= 4 && words[end] == "avec"
                vectors[:,1] .= parse.(Float64, words[1:3])
                hasVector[1] = true
            elseif length(words) >= 4 && words[end] == "bvec"
                vectors[:,2] .= parse.(Float64, words[1:3])
                hasVector[2] = true
            elseif length(words) >= 4 && words[end] == "cvec"
                vectors[:,3] .= parse.(Float64, words[1:3])
                hasVector[3] = true
            elseif length(words) >= 4 && words[end] == "origin"
                origin .= parse.(Float64, words[1:3])
                hasOrigin = true
            elseif words[1] == "Atoms"
                i += 1
                while i <= length(lines)
                    line = strip(lines[i])
                    if isempty(line)
                        i += 1
                        continue
                    end
                    if startswith(line, "#")
                        i += 1
                        continue
                    end
                    words = split(line)
                    atomID = tryparse(Int64, words[1])
                    if atomID === nothing || length(words) < 5
                        break
                    end
                    type = parse(Int64, words[2])
                    x = parse(Float64, words[3])
                    y = parse(Float64, words[4])
                    z = parse(Float64, words[5])
                    push!(types, type)
                    push!(xs, x)
                    push!(ys, y)
                    push!(zs, z)
                    i += 1
                end
                break
            end
            i += 1
        end
    end
    if isempty(types)
        error("No atom coordinates were read from $(fileName).")
    end
    hasOrthBounds = all(isfinite, (xlo, xhi, ylo, yhi, zlo, zhi))
    useVectors = all(hasVector)
    translationVectors = zeros(Float64, 3, 3)
    cellOrigin = zeros(Float64, 3)
    cornerInfoAvailable = false
    if useVectors
        translationVectors .= vectors
        cellOrigin .= hasOrigin ? origin : zeros(Float64, 3)
        cornerInfoAvailable = hasOrigin
    elseif hasOrthBounds
        lengths = [xhi - xlo, yhi - ylo, zhi - zlo]
        if any(l -> l <= 0.0, lengths)
            error("Invalid simulation box bounds read from $(fileName).")
        end
        translationVectors .= [lengths[1] 0.0 0.0; 0.0 lengths[2] 0.0; 0.0 0.0 lengths[3]]
        cellOrigin .= [xlo, ylo, zlo]
        cornerInfoAvailable = true
    else
        error("Could not determine simulation box geometry from $(fileName).")
    end
    rx, ry, rz = replicate
    originalCount = length(types)
    totalAtoms = originalCount * rx * ry * rz
    replicated_types = Vector{Int64}(undef, totalAtoms)
    replicated_xs = Vector{Float64}(undef, totalAtoms)
    replicated_ys = Vector{Float64}(undef, totalAtoms)
    replicated_zs = Vector{Float64}(undef, totalAtoms)
    v1 = translationVectors[:,1]; v2 = translationVectors[:,2]; v3 = translationVectors[:,3]
    idx = 1
    for ix in 0:rx-1
        shift_x_ix = ix * v1[1]
        shift_y_ix = ix * v1[2]
        shift_z_ix = ix * v1[3]
        for iy in 0:ry-1
            shift_x_ixiy = shift_x_ix + iy * v2[1]
            shift_y_ixiy = shift_y_ix + iy * v2[2]
            shift_z_ixiy = shift_z_ix + iy * v2[3]
            for iz in 0:rz-1
                shift_x = shift_x_ixiy + iz * v3[1]
                shift_y = shift_y_ixiy + iz * v3[2]
                shift_z = shift_z_ixiy + iz * v3[3]
                for j in 1:originalCount
                    replicated_types[idx] = types[j]
                    replicated_xs[idx] = xs[j] + shift_x
                    replicated_ys[idx] = ys[j] + shift_y
                    replicated_zs[idx] = zs[j] + shift_z
                    idx += 1
                end
            end
        end
    end
    xmin_atoms = minimum(replicated_xs); xmax_atoms = maximum(replicated_xs)
    ymin_atoms = minimum(replicated_ys); ymax_atoms = maximum(replicated_ys)
    zmin_atoms = minimum(replicated_zs); zmax_atoms = maximum(replicated_zs)
    xmin = xmin_atoms; xmax = xmax_atoms
    ymin = ymin_atoms; ymax = ymax_atoms
    zmin = zmin_atoms; zmax = zmax_atoms
    if cornerInfoAvailable
        total_v1x = v1[1] * rx; total_v1y = v1[2] * rx; total_v1z = v1[3] * rx
        total_v2x = v2[1] * ry; total_v2y = v2[2] * ry; total_v2z = v2[3] * ry
        total_v3x = v3[1] * rz; total_v3y = v3[2] * rz; total_v3z = v3[3] * rz
        corner_xmin = Inf; corner_xmax = -Inf
        corner_ymin = Inf; corner_ymax = -Inf
        corner_zmin = Inf; corner_zmax = -Inf
        for ia in (0, 1)
            for ib in (0, 1)
                for ic in (0, 1)
                    corner_x = cellOrigin[1] + ia * total_v1x + ib * total_v2x + ic * total_v3x
                    corner_y = cellOrigin[2] + ia * total_v1y + ib * total_v2y + ic * total_v3y
                    corner_z = cellOrigin[3] + ia * total_v1z + ib * total_v2z + ic * total_v3z
                    corner_xmin = min(corner_xmin, corner_x)
                    corner_xmax = max(corner_xmax, corner_x)
                    corner_ymin = min(corner_ymin, corner_y)
                    corner_ymax = max(corner_ymax, corner_y)
                    corner_zmin = min(corner_zmin, corner_z)
                    corner_zmax = max(corner_zmax, corner_z)
                end
            end
        end
        xmin = min(xmin, corner_xmin); xmax = max(xmax, corner_xmax)
        ymin = min(ymin, corner_ymin); ymax = max(ymax, corner_ymax)
        zmin = min(zmin, corner_zmin); zmax = max(zmax, corner_zmax)
    end
    replicated_xs .-= xmin
    replicated_ys .-= ymin
    replicated_zs .-= zmin
    xextent = xmax - xmin
    yextent = ymax - ymin
    zextent = zmax - zmin
    if any(extent -> extent <= 0.0, (xextent, yextent, zextent))
        error("Invalid bounding box computed from $(fileName).")
    end
    idx = 1
    return 0.0, xextent, 0.0, yextent, 0.0, zextent, replicated_types, replicated_xs, replicated_ys, replicated_zs
end




function SaveθτData(type_p::Int64, type_t::Int64, 
                       θMatrix::Matrix{Float64}, τMatrix::Matrix{Float64}, 
                       E_p_axis::Vector{Float64}, p_axis::Vector{Float64}, 
                       parameters::Parameters)
    typeDict = parameters.typeDict
    name_p = typeDict[type_p].name
    name_t = typeDict[type_t].name
    fileName = parameters.θτRepository*"/$(name_p)_$(name_t).thetatau"  
    open(fileName, "w") do f
        write(f, "# EPowerRange: $(parameters.EPowerRange)\n")
        write(f, "# pPowerRange: $(parameters.pPowerRange)\n")
        write(f, "# Generated at: $(Dates.now())\n\n")
        write(f, "@ P type: $(name_p) & T type: $(name_t)\n") # to be modified: not neccessry because of the file name has indicated the elements. 
        write(f, "E axis length: $(length(E_p_axis))\n") 
        write(f, "p axis length: $(length(p_axis))\n")
        write(f, "E_p_axis p_axis θ τ\n")
        for (i,E) in enumerate(E_p_axis)
            for (j,p) in enumerate(p_axis)
                write(f,"$(E) $(p) $(θMatrix[i,j]) $(τMatrix[i,j])\n")
            end
        end
    end
end

function parse_range(str)
    parts = split(str, ":")
    start = parse(Float64, parts[1])
    step = parse(Float64, parts[2])
    stop  = parse(Float64, parts[3])
    return start:step:stop
end

function LoadθτData(type_p::Int64, type_t::Int64, parameters::Parameters)
    typeDict = parameters.typeDict
    name_p = typeDict[type_p].name
    name_t = typeDict[type_t].name
    
    E_p_values = Vector{Float64}()
    p_values = Vector{Float64}()
    θ_values = Vector{Float64}()
    τ_values = Vector{Float64}()
    nE = 0
    np = 0
    
    open(parameters.θτRepository*"/$(name_p)_$(name_t).thetatau", "r") do f
        lines = readlines(f)
        i = 1
        while i <= length(lines)
            if startswith(lines[i], "@")
                words = split(lines[i])
                if name_p == words[4] && name_t == words[8]   # to be modified: not neccessry because of the file name has indicated the elements. 
                    i += 1
                    nE = parse(Int64, split(lines[i])[end])
                    i += 1
                    np = parse(Int64, split(lines[i])[end])
                    i += 2
                    while i <= length(lines) && !startswith(lines[i], "@") 
                        if !isempty(lines[i]) && !startswith(lines[i], "#")
                            E_p_value, p_value, θ_value, τ_value = parse.(Float64, split(lines[i]))
                            push!(E_p_values, E_p_value)
                            push!(p_values, p_value)
                            push!(θ_values, θ_value)
                            push!(τ_values, τ_value)
                        end
                        i += 1
                    end              
                    break 
                else
                    i += 1
                end
            elseif startswith(lines[i], "# E_p_power_range")
                word = split(lines[i])[end]
                EPowerRange = parse_range(word)
                if EPowerRange != parameters.EPowerRange
                            log_warning("Warning: The EPowerRange in the file $(name_p)_$(name_t).thetatau is not the same as the EPowerRange in the parameters.")
        log_warning("Generate new thetatau file? (y/n)")
                    answer = readline()
                    if answer == "y"
                        error()
                    end
                end
                i += 1
            elseif startswith(lines[i], "# p_range")
                word = split(lines[i])[end]
                pRange = parse_range(word)
                if pRange != parameters.pRange
                            log_warning("Warning: The pRange in the file $(name_p)_$(name_t).thetatau is not the same as the pRange in the parameters.")
        log_warning("Generate new thetatau file? (y/n)")
                    answer = readline()
                    if answer == "y"
                        error()
                    end
                end
                i += 1
            else
                i += 1
            end
        end
    end
    if length(E_p_values) == 0
        error("Loading θ and τ data for Elements $(name_p) to $(name_t) failed.")
    else
        E_p_axis = sort(unique(E_p_values))
        p_axis = sort(unique(p_values))
        θMatrix = reshape(θ_values, np, nE)' # in julia, data is filled column-wise. 
        τMatrix = reshape(τ_values, np, nE)'
        return E_p_axis, p_axis, θMatrix, τMatrix
    end
end

function LoadDTEData(parameters::Parameters)
    file = parameters.DTEFile
    if !isfile(file)
        error("DTE file $(file) does not exist.")
    end
    DTEData = Vector{Vector{Float64}}()
    environmentCut = 0.0 
    open(file, "r") do f
        lines = readlines(f)
        i = 1
        environmentCut = parse(Float64, split(lines[i])[2])
        while i <= length(lines)
            if startswith(lines[i], "@")
                DTEs = Vector{Float64}()
                i += 1
                while i <= length(lines) && !startswith(lines[i], "@")
                    if !isempty(lines[i]) && !startswith(lines[i], "#")
                        push!(DTEs, parse(Float64,split(lines[i])[1]))
                    end
                    i += 1
                end
                push!(DTEData, DTEs)
            else
                i += 1
            end
        end
    end
    return environmentCut, DTEData
end



function OutputDefects(simulator::Simulator, fileName::String, step::Int64, type::String="a")
    if !simulator.isStore
        error("The defects are not stored, please set isStore to true.")
    end
    open(fileName, type) do file
        write(file, "ITEM: TIMESTEP\n")
        write(file, string(step), "\n")
        write(file, "ITEM: NUMBER OF ATOMS\n")
        write(file, string(length(simulator.displacedAtoms)*2+simulator.numberOfAtoms-simulator.numberOfAtomsWhenStored), "\n")
        write(file, "ITEM: BOX BOUNDS ")
        for d in 1:3
            if simulator.parameters.periodic[d] 
                write(file, "pp ")
            else
                write(file, "ff ")
            end
        end
        write(file, "\n")
        for d in 1:3
            write(file, "0 $(simulator.box.vectors[d,d])\n")
        end
        write(file, "ITEM: ATOMS id type x y z\n")
        for atomIndex in simulator.displacedAtoms
            atom = simulator.atoms[atomIndex]
            write(file, "$(atom.index) $(atom.type) \
            $(atom.coordinate[1]) $(atom.coordinate[2]) $(atom.coordinate[3])\n")
            latticePoint = simulator.latticePoints[atomIndex]  
            write(file, "$(latticePoint.index) $(latticePoint.type+length(keys(simulator.parameters.typeDict))) \
            $(latticePoint.coordinate[1]) $(latticePoint.coordinate[2]) $(latticePoint.coordinate[3])\n")
        end
        for atom in simulator.atoms[simulator.numberOfAtomsWhenStored+1:simulator.numberOfAtoms]
            write(file, "$(atom.index) $(atom.type) \
            $(atom.coordinate[1]) $(atom.coordinate[2]) $(atom.coordinate[3])\n")
        end
    end
end

function OutputAtoms(atoms::Vector{Atom}, simulator::Simulator, fileName::String, step::Int64, type::String="a")
    open(fileName, type) do file
        write(file, "ITEM: TIMESTEP\n")
        write(file, string(step), "\n")
        write(file, "ITEM: NUMBER OF ATOMS\n")
        write(file, string(length(atoms)), "\n")
        write(file, "ITEM: BOX BOUNDS ")
        for d in 1:3
            if simulator.parameters.periodic[d] 
                write(file, "pp ")
            else
                write(file, "ff ")
            end
        end
        write(file, "\n")
        for d in 1:3
            write(file, "0 $(simulator.box.vectors[d,d])\n")
        end
        write(file, "ITEM: ATOMS id type x y z vx vy vz e\n")
        for atom in atoms
            write(file, "$(atom.index) $(atom.type) \
            $(atom.coordinate[1]) $(atom.coordinate[2]) $(atom.coordinate[3]) \
            $(atom.velocityDirection[1]) $(atom.velocityDirection[2]) $(atom.velocityDirection[3]) \
            $(atom.energy)\n")
        end
    end
end

module Output
using Main: Simulator
export @dump, @record

FLUSH_BYTES = 4_096
const _fh = Dict{String, IO}()
const _buf = Dict{String, IOBuffer}()

function _dumpTitle(buf::IOBuffer, simulator::Simulator, step::Int64, atomNumber::Int64)
    print(buf, "ITEM: TIMESTEP\n")
    print(buf, string(step), "\n")
    print(buf, "ITEM: NUMBER OF ATOMS\n")
    print(buf, string(atomNumber), "\n")
    print(buf, "ITEM: BOX BOUNDS ")
    for d in 1:3
        if simulator.parameters.periodic[d] 
            print(buf, "pp ")
        else
            print(buf, "ff ")
        end
    end
    print(buf, "\n")
    for d in 1:3
        print(buf, "0 $(simulator.box.vectors[d,d])\n")
    end
end


function _ensure(file::String, title::String="")
    haskey(_fh, file) && return
    io = open(file, "w")
    _fh[file] = io
    _buf[file] = IOBuffer()
    if title != ""
        print(_buf[file], title, "\n")
    end
end

function _flush!(file::String)
    buf = _buf[file]; io = _fh[file]
    write(io, take!(buf)); flush(io)
end

macro dump(file, atoms, properties=[], stepProperty=:nCascade)
    # for example: properties = ["vx", "vy", "vz", "e"]
    quote
        atomNumber = sum([a.isAlive for a in $(esc(atoms))])
        local _file = $(esc(file))
        Output._ensure(_file)
        local buf = Output._buf[_file]
        Output._dumpTitle(buf, $(esc(:simulator)), $(esc(:simulator)).$(stepProperty), atomNumber)
        print(buf, "ITEM: ATOMS id type x y z ")
        for p in $(esc(properties))
            print(buf, p * " ")
        end
        print(buf, "\n")
        for atom in $(esc(atoms))
            if !atom.isAlive
                continue
            end
            print(buf, string(atom.index) * " " * string(atom.type) * " " * string(atom.coordinate[1]) * " " * string(atom.coordinate[2]) * " " * string(atom.coordinate[3]) * " ")
            for p in $(esc(properties))
                if p[1] == 'v'
                    if p[2] == 'x'
                        print(buf, string(atom.velocityDirection[1]) * " ")
                    elseif p[2] == 'y'
                        print(buf, string(atom.velocityDirection[2]) * " ")
                    elseif p[2] == 'z'
                        print(buf, string(atom.velocityDirection[3]) * " ")
                    else
                        error("Invalid velocity property: $p")
                    end
                elseif p[1] == 'e'
                    print(buf, string(atom.energy) * " ")
                else
                    error("Invalid property: $p")
                end
            end
            print(buf, "\n")
        end
        if position(buf) >= Output.FLUSH_BYTES
            Output._flush!(_file)
        end
    end
end

macro record(file, value, title="", isSmall=false)
    quote
        local _file = $(esc(file))    
        Output._ensure(_file, $(esc(title)))
        local buf = Output._buf[_file]
        print(buf, $(esc(value)), '\n')
        if !$(isSmall)
            if Output.position(buf) ≥ Output.FLUSH_BYTES
                Output._flush!(_file)
            end
        else
            Output._flush!(_file)
        end
    end
end

atexit() do
    for s in keys(_fh)
        _flush!(s); close(_fh[s])
    end
end

end
