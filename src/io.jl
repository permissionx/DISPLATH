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
        write(f, "# pRange: $(parameters.pRange)\n")
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

FLUSH_BYTES = 8_192
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


function _ensure(file::String)
    haskey(_fh, file) && return
    io = open(file, "w")
    _fh[file] = io
    _buf[file] = IOBuffer()
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

macro record(file, value, isSmall=false)
    quote
        local _file = $(esc(file))    
        Output._ensure(_file)
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



# ===== Professional Logging Utilities =====

using Dates

# ANSI color codes
const GRAY = "\033[90m"      # 灰色
const RESET = "\033[0m"      # 重置颜色
const BOLD = "\033[1m"       # 粗体
const DIM = "\033[2m"        # 暗淡

"""
    log_info(args...)

Professional logging with timestamp and gray color formatting.
Replaces standard println with timestamped, colored output.
"""
function log_info(args...)
    timestamp = Dates.format(now(), "HH:MM:SS.sss")
    print(GRAY, "[", timestamp, "] ", RESET)
    
    # Convert all arguments to strings and join them
    message = join(string.(args), "")
    
    # Apply gray color to the message content
    print(GRAY, message, RESET)
    println()
end

"""
    log_success(args...)

Success logging with timestamp and green color.
"""
function log_success(args...)
    timestamp = Dates.format(now(), "HH:MM:SS.sss")
    print("\033[32m", "[", timestamp, "] ✅ ", RESET)  # 绿色
    
    message = join(string.(args), "")
    print("\033[32m", message, RESET)
    println()
end

"""
    log_warning(args...)

Warning logging with timestamp and yellow color.
"""
function log_warning(args...)
    timestamp = Dates.format(now(), "HH:MM:SS.sss")
    print("\033[33m", "[", timestamp, "] ⚠️ ", RESET)  # 黄色
    
    message = join(string.(args), "")
    print("\033[33m", message, RESET)
    println()
end

"""
    log_error(args...)

Error logging with timestamp and red color.
"""
function log_error(args...)
    timestamp = Dates.format(now(), "HH:MM:SS.sss")
    print("\033[31m", "[", timestamp, "] ❌ ", RESET)  # 红色
    
    message = join(string.(args), "")
    print("\033[31m", message, RESET)
    println()
end

"""
    log_debug(args...)

Debug logging with timestamp and dim gray color.
"""
function log_debug(args...)
    timestamp = Dates.format(now(), "HH:MM:SS.sss")
    print(DIM, GRAY, "[", timestamp, "] 🐛 ", RESET)
    
    message = join(string.(args), "")
    print(DIM, GRAY, message, RESET)
    println()
end

# Export the logging functions
export log_info, log_success, log_warning, log_error, log_debug
