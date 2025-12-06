#将模拟器状态输出为LAMMPS格式
function Dump(simulator::Simulator, fileName::String, step::Int64, type::String="a", isDebug::Bool=false)# ... 原子数据输出实现
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
                    write(
                        file,
                        "$(atom.index) $(atom.type) \
            $(atom.coordinate[1]) $(atom.coordinate[2]) $(atom.coordinate[3]) \
            $(atom.velocityDirection[1]*sqrt(2*atom.mass*atom.energy)) $(atom.velocityDirection[2]*sqrt(2*atom.mass*atom.energy)) $(atom.velocityDirection[3]*sqrt(2*atom.mass*atom.energy)) \
            $(atom.energy) \
            $(atom.cellIndex[1]) $(atom.cellIndex[2]) $(atom.cellIndex[3]) \
            $(GetDTE(atom, simulator))\n"
                    )
                else
                    write(
                        file,
                        "$(atom.index) $(atom.type) \
            $(atom.coordinate[1]) $(atom.coordinate[2]) $(atom.coordinate[3]) $(atom.energy)\n"
                    )
                end
            end
        end
    end
end

function ReadDate(fileName::String, replicate::Vector{Int64})  # 修改1: 新增replicate参数，支持晶格复制
    """
    从数据文件读取原子信息，支持非正交盒子（通过读取基矢）和晶格复制功能。
    
    参数:
    - fileName: 输入文件名
    - replicate: 三维复制因子，例如[2,2,1]表示在x,y方向各复制2次，z方向不复制
    
    返回:
    - xlo, xhi, ylo, yhi, zlo, zhi: 归一化后的盒子边界（原子坐标平移至原点开始）
    - types: 原子类型数组
    - xs, ys, zs: 归一化后的原子坐标数组
    
    功能增强:
    1. 支持读取非正交盒子的基矢 (avec, bvec, cvec) 和原点 (origin)
    2. 支持通过复制因子扩大模拟体系
    3. 更健壮的错误处理和边界情况检测
    """
    if length(replicate) != 3
        error("Replicate vector must contain 3 entries for x/y/z directions.")
    end
    if any(r -> r < 1, replicate)
        error("Replicate factors must be >= 1.")
    end
    
    # 初始化变量 - 使用NaN表示未初始化状态，提高代码健壮性
    xlo = NaN; xhi = NaN
    ylo = NaN; yhi = NaN
    zlo = NaN; zhi = NaN
    vectors = zeros(Float64, 3, 3)  # 存储盒子基矢矩阵 (3×3)
    hasVector = falses(3)           # 标记基矢是否被读取
    origin = zeros(Float64, 3)      # 盒子原点坐标
    hasOrigin = false               # 标记原点是否被读取
    
    types = Int64[]
    xs = Float64[]
    ys = Float64[]
    zs = Float64[]
    
    open(fileName, "r") do f
        lines = readlines(f)
        i = 1
        while i <= length(lines)
            line = strip(lines[i])  # 去除首尾空白字符，提高解析鲁棒性
            if isempty(line) || startswith(line, "#")
                i += 1
                continue
            end
            
            words = split(line)
            
            # 读取正交盒子边界（传统LAMMPS格式）
            if length(words) >= 4 && words[4] == "xhi"
                xlo = parse(Float64, words[1])
                xhi = parse(Float64, words[2])
            elseif length(words) >= 4 && words[4] == "yhi"
                ylo = parse(Float64, words[1])
                yhi = parse(Float64, words[2])
            elseif length(words) >= 4 && words[4] == "zhi"
                zlo = parse(Float64, words[1])
                zhi = parse(Float64, words[2])
            
            # 读取非正交盒子基矢（扩展格式支持）
            elseif length(words) >= 4 && words[end] == "avec"
                vectors[:,1] .= parse.(Float64, words[1:3])
                hasVector[1] = true
            elseif length(words) >= 4 && words[end] == "bvec"
                vectors[:,2] .= parse.(Float64, words[1:3])
                hasVector[2] = true
            elseif length(words) >= 4 && words[end] == "cvec"
                vectors[:,3] .= parse.(Float64, words[1:3])
                hasVector[3] = true
            
            # 读取盒子原点坐标
            elseif length(words) >= 4 && words[end] == "origin"
                origin .= parse.(Float64, words[1:3])
                hasOrigin = true
            
            # 读取原子数据部分
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
                    # 原子ID解析失败或数据不足时退出原子读取循环
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
                break  # 原子数据读取完成后退出主循环
            end
            i += 1
        end
    end
    
    # 验证是否成功读取原子数据
    if isempty(types)
        error("No atom coordinates were read from $(fileName).")
    end
    
    # 确定盒子几何信息的来源：优先使用基矢，其次使用正交边界
    hasOrthBounds = all(isfinite, (xlo, xhi, ylo, yhi, zlo, zhi))
    useVectors = all(hasVector)  # 三个基矢都齐全时使用基矢表示
    
    translationVectors = zeros(Float64, 3, 3)
    cellOrigin = zeros(Float64, 3)
    cornerInfoAvailable = false  # 标记是否有完整的盒子角点信息
    
    if useVectors
        # 使用基矢定义的非正交盒子
        translationVectors .= vectors
        cellOrigin .= hasOrigin ? origin : zeros(Float64, 3)
        cornerInfoAvailable = hasOrigin
    elseif hasOrthBounds
        # 使用传统正交盒子边界
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
    
    # 执行晶格复制操作
    rx, ry, rz = replicate
    originalCount = length(types)
    totalAtoms = originalCount * rx * ry * rz
    
    # 预分配数组以提高性能
    replicated_types = Vector{Int64}(undef, totalAtoms)
    replicated_xs = Vector{Float64}(undef, totalAtoms)
    replicated_ys = Vector{Float64}(undef, totalAtoms)
    replicated_zs = Vector{Float64}(undef, totalAtoms)
    
    v1 = translationVectors[:,1]
    v2 = translationVectors[:,2]
    v3 = translationVectors[:,3]
    
    idx = 1
    for ix in 0:rx-1
        # 计算x方向复制引起的位移
        shift_x_ix = ix * v1[1]
        shift_y_ix = ix * v1[2]
        shift_z_ix = ix * v1[3]
        
        for iy in 0:ry-1
            # 叠加y方向复制引起的位移
            shift_x_ixiy = shift_x_ix + iy * v2[1]
            shift_y_ixiy = shift_y_ix + iy * v2[2]
            shift_z_ixiy = shift_z_ix + iy * v2[3]
            
            for iz in 0:rz-1
                # 叠加z方向复制引起的位移
                shift_x = shift_x_ixiy + iz * v3[1]
                shift_y = shift_y_ixiy + iz * v3[2]
                shift_z = shift_z_ixiy + iz * v3[3]
                
                # 复制原始晶格中的所有原子
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
    
    # 计算复制后所有原子的边界
    xmin_atoms = minimum(replicated_xs); xmax_atoms = maximum(replicated_xs)
    ymin_atoms = minimum(replicated_ys); ymax_atoms = maximum(replicated_ys)
    zmin_atoms = minimum(replicated_zs); zmax_atoms = maximum(replicated_zs)
    
    xmin = xmin_atoms; xmax = xmax_atoms
    ymin = ymin_atoms; ymax = ymax_atoms
    zmin = zmin_atoms; zmax = zmax_atoms
    
    # 如果有盒子角点信息，计算完整的体系边界（包括原子和盒子边界）
    if cornerInfoAvailable
        # 计算复制后的总基矢
        total_v1x = v1[1] * rx; total_v1y = v1[2] * rx; total_v1z = v1[3] * rx
        total_v2x = v2[1] * ry; total_v2y = v2[2] * ry; total_v2z = v2[3] * ry
        total_v3x = v3[1] * rz; total_v3y = v3[2] * rz; total_v3z = v3[3] * rz
        
        # 遍历所有8个角点，找到最小和最大坐标
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
        
        # 取原子边界和盒子边界的并集
        xmin = min(xmin, corner_xmin); xmax = max(xmax, corner_xmax)
        ymin = min(ymin, corner_ymin); ymax = max(ymax, corner_ymax)
        zmin = min(zmin, corner_zmin); zmax = max(zmax, corner_zmax)
    end
    
    # 归一化坐标：将整个体系平移，使最小坐标为0
    replicated_xs .-= xmin
    replicated_ys .-= ymin
    replicated_zs .-= zmin
    
    # 计算归一化后的盒子尺寸
    xextent = xmax - xmin
    yextent = ymax - ymin
    zextent = zmax - zmin
    
    # 验证边界有效性
    if any(extent -> extent <= 0.0, (xextent, yextent, zextent))
        error("Invalid bounding box computed from $(fileName).")
    end
    
    # 返回归一化后的边界和坐标
    return 0.0, xextent, 0.0, yextent, 0.0, zextent, replicated_types, replicated_xs, replicated_ys, replicated_zs
end



#θτ数据的保存和加载
function SaveθτData(type_p::Int64, type_t::Int64,
    θMatrix::Matrix{Float64}, τMatrix::Matrix{Float64},
    E_p_axis::Vector{Float64}, p_axis::Vector{Float64},
    parameters::Parameters)
    typeDict = parameters.typeDict
    name_p = typeDict[type_p].name
    name_t = typeDict[type_t].name
    fileName = parameters.θτRepository * "/$(name_p)_$(name_t).thetatau"
    open(fileName, "w") do f
        write(f, "# EPowerRange: $(parameters.EPowerRange)\n")
        write(f, "# pPowerRange: $(parameters.pPowerRange)\n")
        write(f, "# Generated at: $(Dates.now())\n\n")
        write(f, "@ P type: $(name_p) & T type: $(name_t)\n") # to be modified: not neccessry because of the file name has indicated the elements. 
        write(f, "E axis length: $(length(E_p_axis))\n")
        write(f, "p axis length: $(length(p_axis))\n")
        write(f, "E_p_axis p_axis θ τ\n")
        for (i, E) in enumerate(E_p_axis)
            for (j, p) in enumerate(p_axis)
                write(f, "$(E) $(p) $(θMatrix[i,j]) $(τMatrix[i,j])\n")
            end
        end
    end
end

function parse_range(str)
    parts = split(str, ":")
    start = parse(Float64, parts[1])
    step = parse(Float64, parts[2])
    stop = parse(Float64, parts[3])
    return start:step:stop
end

#θτ数据的保存和加载
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

    open(parameters.θτRepository * "/$(name_p)_$(name_t).thetatau", "r") do f #存储预先计算的散射角(θ)和飞行时间(τ)数据
        lines = readlines(f)
        i = 1
        while i <= length(lines)
            if startswith(lines[i], "@")
                words = split(lines[i])
                if name_p == words[4] && name_t == words[8]   # to be modified: not neccessry because of the file name has indicated the elements. 
                    i += 1
                    nE = parse(Int64, split(lines[i])[end])# 读取能量轴长度
                    i += 1
                    np = parse(Int64, split(lines[i])[end])# 读取碰撞参数轴长度
                    i += 2  # 跳过表头
                    while i <= length(lines) && !startswith(lines[i], "@") #数据读取循环
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
            elseif startswith(lines[i], "# E_p_power_range") #能量范围验证
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

function LoadDTEData(parameters::Parameters)# 动态时间演化数据加载
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
                        push!(DTEs, parse(Float64, split(lines[i])[1]))
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



function OutputDefects(simulator::Simulator, fileName::String, step::Int64, type::String="a")#缺陷统计输出
    if !simulator.isStore
        error("The defects are not stored, please set isStore to true.")
    end
    open(fileName, type) do file
        write(file, "ITEM: TIMESTEP\n")
        write(file, string(step), "\n")
        write(file, "ITEM: NUMBER OF ATOMS\n")
        write(file, string(length(simulator.displacedAtoms) * 2 + simulator.numberOfAtoms - simulator.numberOfAtomsWhenStored), "\n")
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
            write(
                file,
                "$(atom.index) $(atom.type) \
    $(atom.coordinate[1]) $(atom.coordinate[2]) $(atom.coordinate[3])\n"
            )
            latticePoint = simulator.latticePoints[atomIndex]
            write(
                file,
                "$(latticePoint.index) $(latticePoint.type+length(keys(simulator.parameters.typeDict))) \
    $(latticePoint.coordinate[1]) $(latticePoint.coordinate[2]) $(latticePoint.coordinate[3])\n"
            )
        end
        for atom in simulator.atoms[simulator.numberOfAtomsWhenStored+1:simulator.numberOfAtoms]
            write(
                file,
                "$(atom.index) $(atom.type) \
    $(atom.coordinate[1]) $(atom.coordinate[2]) $(atom.coordinate[3])\n"
            )
        end
    end
end

function OutputAtoms(atoms::Vector{Atom}, simulator::Simulator, fileName::String, step::Int64, type::String="a")#原子数据输出
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
            write(
                file,
                "$(atom.index) $(atom.type) \
    $(atom.coordinate[1]) $(atom.coordinate[2]) $(atom.coordinate[3]) \
    $(atom.velocityDirection[1]) $(atom.velocityDirection[2]) $(atom.velocityDirection[3]) \
    $(atom.energy)\n"
            )
        end
    end
end

module Output
using Main: Simulator
export @dump, @record

FLUSH_BYTES = 4_096
const _fh = Dict{String,IO}()
const _buf = Dict{String,IOBuffer}()

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


function _ensure(file::String, title::String="") #懒初始化文件句柄
    haskey(_fh, file) && return  # 文件已初始化则直接返回
    io = open(file, "w")          # 创建文件
    _fh[file] = io               # 存储文件句柄
    _buf[file] = IOBuffer()      # 创建缓冲区
    if title != ""
        print(_buf[file], title, "\n")  # 写入表头
    end
end

function _flush!(file::String)
    buf = _buf[file]; io = _fh[file]
    write(io, take!(buf))  # 将缓冲区内容写入文件
    flush(io)              # 强制刷新操作系统缓冲区
end
#file: 输出文件名   atoms: 原子集合（Vector{Atom}）  properties: 可选属性列表，如 ["vx", "vy", "vz", "e"]  stepProperty: 时间步属性，默认使用simulator.nCascade
# for example: properties = ["vx", "vy", "vz", "e"]
macro dump(file, atoms, properties=[], stepProperty=:nCascade)
    quote
        # 计算当前存活的原子总数
        # 通过列表推导式遍历所有原子，检查每个原子的isAlive状态
        # sum函数统计所有存活原子（isAlive=true）的数量
        # 这个计数用于确定dump文件中要写入的原子数量
        atomNumber = sum([a.isAlive for a in $(esc(atoms))])

        # 将文件名字符串转换为局部变量，确保宏卫生性
        # esc()函数用于转义符号，防止宏展开时的变量捕获问题
        # _file变量存储输出文件的路径和名称
        local _file = $(esc(file))

        # 懒初始化文件句柄和缓冲区
        # 如果该文件尚未在全局字典中注册，则创建对应的文件句柄和IO缓冲区
        # 这个机制确保文件只在第一次使用时打开，提高性能
        Output._ensure(_file)

        # 获取对应文件的IO缓冲区引用
        # 使用文件名字符串作为键，从全局缓冲区字典中获取对应的IOBuffer
        # 所有原子数据首先写入这个内存缓冲区，积累到一定量后再刷新到磁盘
        local buf = Output._buf[_file]

        # 写入LAMMPS兼容的格式头信息到缓冲区
        # 这个函数生成标准的LAMMPS dump文件头部，包括：
        # - TIMESTEP: 当前时间步长
        # - NUMBER OF ATOMS: 原子总数
        # - BOX BOUNDS: 模拟盒子的边界条件和尺寸
        # - ATOMS: 原子数据的列标题
        # 使用simulator对象获取盒子信息和边界条件
        # stepProperty指定从simulator的哪个字段获取当前时间步（如nCascade）
        Output._dumpTitle(buf, $(esc(:simulator)), $(esc(:simulator)).$(stepProperty), atomNumber)
        print(buf, "ITEM: ATOMS id type x y z ")#列标题生成
        for p in $(esc(properties))
            print(buf, p * " ")
        end
        print(buf, "\n")
        for atom in $(esc(atoms)) #原子数据输出
            if !atom.isAlive
                continue
            end
            # 输出基础属性
            print(buf, string(atom.index) * " " * string(atom.type) * " " * string(atom.coordinate[1]) * " " * string(atom.coordinate[2]) * " " * string(atom.coordinate[3]) * " ")
            # 输出额外属性
            for p in $(esc(properties))
                if p[1] == 'v' # 速度属性
                    if p[2] == 'x'
                        print(buf, string(atom.velocityDirection[1]) * " ")
                    elseif p[2] == 'y'
                        print(buf, string(atom.velocityDirection[2]) * " ")
                    elseif p[2] == 'z'
                        print(buf, string(atom.velocityDirection[3]) * " ")
                    else
                        error("Invalid velocity property: $p")
                    end
                elseif p[1] == 'e' # 能量属性
                    print(buf, string(atom.energy) * " ")
                else
                    error("Invalid property: $p")
                end
            end
            print(buf, "\n")
        end
        if position(buf) >= Output.FLUSH_BYTES #缓冲刷新
            Output._flush!(_file)
        end
    end
end

macro record(file, value, title="", isSmall=false)#@record宏应用：CSV数据记录 日志输出 小文件立即刷新选项
    #宏参数：file：输出文件名;value：要记录的数据值;title：可选的标题/表头（默认为空）;isSmall：刷新策略标志（默认为false）
    quote
        local _file = $(esc(file))#创建局部变量_file存储文件名,$(esc(file))：对file参数进行卫生处理，避免变量捕获问题,esc确保宏内部的变量不会与外部的同名变量冲突
        Output._ensure(_file, $(esc(title)))#调用Output._ensure函数确保文件存在且已正确初始化;如果文件不存在，创建文件并写入title作为表头;如果文件已存在，检查是否已有相同表头
        local buf = Output._buf[_file] #作用：从全局字典Output._buf中获取对应文件的缓冲区;每个文件有独立的IO缓冲区，避免多文件写入冲突
        print(buf, $(esc(value)), '\n')#操作：将转义后的value和换行符写入缓冲区;使用print而不是println因为显式添加了'\n';缓冲区在内存中累积数据，减少磁盘I/O操作
        if !$(isSmall) #刷新策略控制-双重刷新策略：A.大文件模式 (isSmall=false)
            #机制：检查缓冲区当前位置是否达到刷新阈值FLUSH_BYTES（4096字节）;达到阈值时才执行实际的文件写入;优点：减少磁盘I/O，提高大批量数据写入性能
            if Output.position(buf) ≥ Output.FLUSH_BYTES
                Output._flush!(_file)
            end
        else#B. 小文件模式 (isSmall=true).机制：立即刷新缓冲区，确保数据实时写入磁盘;优点：数据不会因程序崩溃而丢失，适合关键日志
            Output._flush!(_file)
        end
    end
end

#优雅退出：程序退出时自动刷新所有缓冲区 关闭所有文件句柄 防止数据丢失
atexit() do
    for s in keys(_fh)
        _flush!(s)
        close(_fh[s])
    end
end

end

