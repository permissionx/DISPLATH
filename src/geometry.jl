# 导入StaticArrays包，用于高性能的静态数组操作
using StaticArrays

# 定义Box构造函数：根据输入向量创建模拟盒子
function Box(Vectors::Matrix{Float64})
    # 输入验证
    if size(Vectors) != (3, 3)
        throw(GeometryError("Box construction", "Vectors must be a 3×3 matrix, got $(size(Vectors))"))
    end
    det_vectors = det(Vectors)
    if det_vectors <= 0.0
        throw(GeometryError("Box construction", "Vectors determinant must be positive (volume > 0), got $(det_vectors)"))
    end
    
    # 记录盒子创建信息，显示盒子的三维尺寸（四舍五入到2位小数）
    log_info("Box created: $(round(Vectors[1,1]; digits=2)) × $(round(Vectors[2,2]; digits=2)) × $(round(Vectors[3,3]; digits=2)) Å")
    # 返回Box结构体：包含原始向量、逆矩阵转置（倒易空间向量）、正交标志设为true
    return Box(Vectors, inv(Vectors'), true)
end 

# 根据原胞基向量和尺寸创建模拟盒子
function CreateBoxByPrimaryVectors(primaryVectors::Matrix{Float64}, sizes::Vector{Int64})
    # 将原胞基向量按尺寸放大得到盒子基向量（单位：Å）
    vectors = primaryVectors .* sizes
    # 调用Box构造函数创建最终的模拟盒子
    return Box(vectors)
end 

# 原子构造函数：根据类型、坐标和参数创建原子对象
function Atom(type::Int64, coordinate::Vector{Float64}, parameters::Parameters)
    # 输入验证
    if type <= 0
        throw(InvalidParameterError("type", type, "Atom type must be a positive integer"))
    end
    if !haskey(parameters.typeDict, type)
        throw(InvalidParameterError("type", type, "Atom type not found in typeDict"))
    end
    if length(coordinate) != 3
        throw(GeometryError("Atom construction", "Coordinate must be a 3-element vector, got length $(length(coordinate))"))
    end
    if any(!isfinite, coordinate)
        throw(GeometryError("Atom construction", "Coordinate contains non-finite values: $(coordinate)"))
    end
    
    # 初始化原子基本属性
    index = 0                         # 原子索引，初始为0，后续会分配唯一值
    isAlive = true                    # 原子存活状态，初始为true
    cellIndex = (0,0,0)               # 所在网格单元索引，初始化为(0,0,0)
    velocityDirection = SVector{3,Float64}(0.0, 0.0, 0.0)  # 速度方向向量，初始为零向量
    energy = 0.0                      # 原子能量，初始为0 eV
    
    # 从类型字典中获取原子的物理属性
    radius, mass, Z, dte, bde, _, _ = TypeToProperties(type, parameters.typeDict)
    
    # 碰撞相关参数初始化
    #numberOfEmptyCells = 0           # 注释：穿越的空单元数量
    emptyPath = 0.0                   # 自由飞行路径长度，初始为0 Å
    pValue = 0.0                      # 碰撞参数值，初始为0 Å
    pVector = SVector{3,Float64}(0.0, 0.0, 0.0)   # 碰撞向量，初始为零向量
    pPoint = SVector{3,Float64}(0.0, 0.0, 0.0)    # 碰撞点坐标，初始为零向量
    lastTargets = Vector{Int64}()     # 上一次碰撞的目标原子索引列表，初始为空
    pL = 0.0                          # 碰撞路径长度参数，初始为0 Å
    pAtomIndex = -1                   # 临时字段：关联的入射原子索引，初始为-1
    pDirection = Float64[0.0,0.0,0.0] # 临时字段：碰撞方向向量，初始为零向量
    
    # 晶格相关参数初始化
    latticePointIndex = -1             # 晶格点索引，-1表示不在晶格位置上
    
    # KMC模拟相关参数初始化
    frequency = 0.0                   # 迁移频率，初始为0 Hz
    frequencies = Vector{Float64}()   # 各迁移路径的频率列表，初始为空
    finalLatticePointEnvIndexs = Vector{Int64}()  # 目标晶格点环境索引列表，初始为空
    eventIndex = -1                   # KMC事件索引，初始为-1
    
    # 动态加载相关参数初始化
    isNewlyLoaded = false             # 新加载标志，初始为false
    lattcieCoordinate = SVector{3,Float64}(coordinate[1], coordinate[2], coordinate[3])  # 理想晶格坐标
    indexInCell = 0                   # 在单元内的索引，初始为0
    
    # 创建并返回完整的Atom结构体
    return Atom(index, isAlive, type, coordinate[:], cellIndex, 
                radius, mass, velocityDirection, energy, Z, 
                dte, bde, emptyPath, #numberOfEmptyCells,
                pValue, pVector, pPoint, pL, pAtomIndex, pDirection, lastTargets, # temperory 
                latticePointIndex,
                frequency, frequencies, finalLatticePointEnvIndexs, eventIndex, 
                isNewlyLoaded, lattcieCoordinate, indexInCell)
end

# 根据原子类型从类型字典中获取物理属性
function TypeToProperties(type::Int64, typeDict::Dict{Int64, Element})
    # 检查类型是否存在于字典中
    if haskey(typeDict, type)
        # 获取对应的元素对象
        element = typeDict[type]
        # 返回元素的物理属性：半径、质量、原子序数、位移阈值能量、结合能、α参数、β参数
        return element.radius, element.mass, element.Z, element.dte, element.bde, element.alpha, element.beta 
    else
        # 如果类型不存在，抛出错误
        throw(InvalidParameterError("type", type, "Atom type not found in typeDict"))
    end 
end 

# 设置网格单元的邻近单元信息（优化版本，使用三重循环而非递归）
function SetCellNeighborInfo!(cell::Cell, grid::Grid)
    # 遍历所有可能的邻近单元偏移量（-1, 0, 1）在每个维度
    for delta_x in [-1, 0, 1]
        for delta_y in [-1, 0, 1]
            for delta_z in [-1, 0, 1]
                # 将偏移量转换为Int8元组
                neighborKeys = (Int8(delta_x), Int8(delta_y), Int8(delta_z))  
                neighborIndex = [0, 0, 0]  # 初始化邻近单元索引
                neighborCross = [Int8(0), Int8(0), Int8(0)]  # 初始化边界穿越标志
                
                # 对每个维度计算邻近单元索引和边界穿越标志
                for d in 1:3
                    delta = neighborKeys[d]  # 当前维度的偏移量
                    index = cell.index[d] + delta  # 计算原始索引
                    cross = Int8(0)  # 初始化穿越标志
                    
                    # 处理周期性边界条件：索引越界时的回绕
                    if index < 1
                        index += grid.sizes[d]  # 负向越界，正向回绕
                        cross = Int8(-1)        # 设置负向穿越标志
                    elseif index > grid.sizes[d]
                        index -= grid.sizes[d]  # 正向越界，负向回绕  
                        cross = Int8(1)         # 设置正向穿越标志
                    end
                    neighborIndex[d] = index      # 存储计算后的索引
                    neighborCross[d] = cross      # 存储穿越标志
                end
                
                # 将索引和穿越标志转换为元组格式
                neighborIndex_tuple = (neighborIndex[1], neighborIndex[2], neighborIndex[3])
                neighborCross_tuple = (neighborCross[1], neighborCross[2], neighborCross[3])
                
                # 创建邻近单元信息结构体
                neighborCellInfo = NeighborCellInfo(neighborIndex_tuple, neighborCross_tuple)
                
                # 计算在neighborCellsInfo数组中的索引（将偏移量从[-1,0,1]映射到[1,2,3]）
                idx = (delta_x + 2, delta_y + 2, delta_z + 2)
                
                # 将邻近单元信息存储到对应的数组位置
                cell.neighborCellsInfo[idx...] = neighborCellInfo
            end
        end
    end
end

# 根据盒子和输入向量创建网格系统
function CreateGrid(box::Box, inputVectors::Matrix{Float64}, parameters::Parameters)
    # 检查盒子是否为正交系统，非正交系统暂不支持
    if !box.isOrthogonal
        throw(GeometryError("CreateGrid", "The box is not orthogonal, please use the orthogonal box."))
    end
    
    # 初始化网格尺寸和向量数组
    sizes = Vector{Int64}(undef, 3)    # 网格各维度的单元数量
    vectors = Matrix{Float64}(undef, 3, 3)  # 单个网格单元的基向量
    
    # 计算每个维度的网格参数
    for d in 1:3
        # 计算该维度的网格单元数量（取整）
        sizes[d] = Int64(floor(box.vectors[d,d] / inputVectors[d,d]))
        # 计算该维度网格单元的实际大小
        vectors[d,d] = box.vectors[d,d] / sizes[d]
    end
    
    # 记录网格创建信息
    log_info("Cell grid: $(sizes[1]) × $(sizes[2]) × $(sizes[3]) = $(sizes[1]*sizes[2]*sizes[3]) cells")
    log_info("Cell size: $(round(vectors[1,1]; digits=2)) × $(round(vectors[2,2]; digits=2)) × $(round(vectors[3,3]; digits=2)) Å")
    
    # 根据动态加载标志选择不同的网格存储方式
    if !parameters.is_dynamic_load
        # 静态加载：使用三维数组存储网格单元
        cells = Array{Cell, 3}(undef, sizes[1], sizes[2], sizes[3])
        # 显示进度条创建所有网格单元
        @showprogress desc="Creating cells: " for x in 1:sizes[1]
            for y in 1:sizes[2]
                for z in 1:sizes[3]
                    # 为每个网格位置创建单元
                    cells[x, y, z] = CreateCell((x, y, z), vectors)
                end
            end    
        end
        # 计算单个网格单元的体积
        cellVolume = vectors[1,1] * vectors[2,2] * vectors[3,3]
        # 创建网格结构体
        grid = Grid(cells, vectors, sizes, cellVolume) 
        # 显示进度条初始化所有单元的邻近信息
        @showprogress desc="Pushing cell neighbors: " for cell in grid.cells
            SetCellNeighborInfo!(cell, grid)
        end
    else
        # 动态加载：使用字典存储网格单元，支持稀疏网格
        cells = Dict{Tuple{Int64, Int64, Int64}, Cell}()    
        # 计算单个网格单元的体积
        cellVolume = vectors[1,1] * vectors[2,2] * vectors[3,3]
        # 创建网格结构体（邻近信息在需要时动态初始化）
        grid = Grid(cells, vectors, sizes, cellVolume) 
    end
    
    # 记录网格创建成功信息
    log_success("Cell grid created")
    log_separator()
    return grid
end

# 静态加载模式下获取网格单元（三维数组访问）
function _GetCellDense(grid::Grid, cellIndex::Tuple{Int64, Int64, Int64})
    x, y, z = cellIndex
    # 边界检查
    if x < 1 || x > grid.sizes[1] || y < 1 || y > grid.sizes[2] || z < 1 || z > grid.sizes[3]
        throw(GeometryError("GetCell", "Cell index $cellIndex is out of bounds. Grid size: $(grid.sizes)"))
    end
    return grid.cells[x, y, z]  # 使用索引访问三维数组
end

# 创建单个网格单元
function CreateCell(cellIndex::Tuple{Int64, Int64, Int64}, vectors::Matrix{Float64})
    x, y, z = cellIndex  # 解包网格索引
    
    # 计算网格单元的空间边界范围
    ranges = Matrix{Float64}(undef, 3, 2)
    ranges[1,1] = (x-1) * vectors[1,1]  # x方向下边界
    ranges[1,2] = x * vectors[1,1]      # x方向上边界  
    ranges[2,1] = (y-1) * vectors[2,2]  # y方向下边界
    ranges[2,2] = y * vectors[2,2]      # y方向上边界
    ranges[3,1] = (z-1) * vectors[3,3]  # z方向下边界
    ranges[3,2] = z * vectors[3,3]      # z方向上边界
    
    # 创建并返回网格单元结构体，初始化所有字段
    cell = Cell(cellIndex, Vector{Atom}(), Vector{LatticePoint}(), 
                            ranges, 
                            Array{NeighborCellInfo, 3}(undef, 3, 3, 3), false, 0.0)
    return cell
end

# 动态加载模式下获取或创建网格单元（字典访问）
function _GetCellDict!(grid::Grid, cellIndex::Tuple{Int64, Int64, Int64})
    # 使用get!函数：如果键存在则返回值，否则创建新值
    return get!(grid.cells, cellIndex) do
        CreateCell(cellIndex, grid.vectors)  # 创建新单元的闭包
    end
end 

# 统一的网格单元获取接口，根据加载模式选择实现
# 注意：此函数需要根据 grid.cells 的实际类型自动判断
function GetCell(grid::Grid, cellIndex::Tuple{Int64, Int64, Int64})
    # 基本输入验证
    x, y, z = cellIndex
    if x < 1 || y < 1 || z < 1
        throw(GeometryError("GetCell", "Cell index components must be positive, got $cellIndex"))
    end
    
    if isa(grid.cells, Array)
        return _GetCellDense(grid, cellIndex)  # 静态加载：数组访问（内部有边界检查）
    else
        return _GetCellDict!(grid, cellIndex)  # 动态加载：字典访问
    end
end

# 初始化类型相关的物理常数
function InitConstantsByType(typeDict::Dict{Int64, Element}, parameters::Parameters)
    # 初始化各种常数字典
    V_upterm = Dict{Vector{Int64}, Float64}()      # 势函数上项常数
    a_U = Dict{Vector{Int64}, Float64}()           # 通用屏蔽长度
    E_m = Dict{Int64, Float64}()                   # 电子阻止特征能量
    S_e_upTerm = Dict{Vector{Int64}, Float64}()    # 电子阻止上项
    S_e_downTerm = Dict{Vector{Int64}, Float64}()  # 电子阻止下项
    x_nl = Dict{Vector{Int64}, Float64}()          # 非局域距离参数
    a = Dict{Vector{Int64}, Float64}()             # 指数衰减长度
    Q_nl = Dict{Vector{Int64}, Float64}()          # 非局域能量损失系数
    Q_loc = Dict{Vector{Int64}, Float64}()         # 局域能量损失系数
    types = keys(typeDict)                         # 所有原子类型
    qMax = Dict{Vector{Int64}, Float64}()          # 最大碰撞参数
    sigma = Dict{Int64, Float64}()                 # 热振动均方根位移
    
    # 记录振动参数信息
    log_info("")
    log_info("Vibration σ for each type:")
    
    # 遍历所有原子类型对，计算类型相关的常数
    for p in types
        # 获取入射原子类型的物理属性
        radius_p, mass_p, Z_p, _, _, α_p, β_p = TypeToProperties(p, typeDict)
        
        # 对每个目标原子类型计算碰撞相关的常数
        for t in types
            # 获取目标原子类型的物理属性（只需要半径和原子序数）
            radius_t, _, Z_t, _, _, _, _ = TypeToProperties(t, typeDict)
            
            # 计算各种物理常数（通过BCA模块的常数函数）
            V_upterm[[p,t]] = BCA.ConstantFunctions.V_upterm(Z_p, Z_t)                    # 势函数上项
            a_U[[p,t]] = BCA.ConstantFunctions.a_U(Z_p, Z_t)                             # 通用屏蔽长度
            S_e_upTerm[[p,t]] = BCA.ConstantFunctions.S_e_upTerm(p, Z_p, Z_t, mass_p, α_p) # 电子阻止上项
            x_nl[[p,t]] = BCA.ConstantFunctions.x_nl(p, Z_p, Z_t, β_p)                    # 非局域距离参数
            a[[p,t]] = BCA.ConstantFunctions.a(Z_p, Z_t)                                 # 指数衰减长度
            Q_nl[[p,t]] = BCA.ConstantFunctions.Q_nl(Z_p, Z_t, parameters.pMax)          # 非局域能量损失系数
            Q_loc[[p,t]] = BCA.ConstantFunctions.Q_loc(Z_p, Z_t)                         # 局域能量损失系数
            qMax[[p,t]] = radius_p + radius_t                                            # 最大碰撞参数（原子半径和）
        end
        
        # 计算单个原子类型的常数
        E_m[p] = BCA.ConstantFunctions.E_m(Z_p, mass_p)  # 电子阻止特征能量
        sigma[p] = TemperatureToSigma(parameters.temperature, parameters.DebyeTemperature, mass_p)  # 热振动参数
        
        # 记录该原子类型的振动参数
        log_info("  Type $(p): σ = $(round(sigma[p]; digits=3)) Å")
    end
    
    # 返回完整的常数集合结构体
    return ConstantsByType(V_upterm, a_U, E_m, S_e_upTerm, S_e_downTerm, x_nl, a, Q_nl, Q_loc, qMax, sigma)
end

# 初始化θ和τ插值函数，用于快速计算碰撞过程中的散射角和飞行时间参数
function InitθτFunctions(parameters::Parameters, constantsByType::ConstantsByType)
    typeDict = parameters.typeDict  # 获取原子类型字典
    θFunctions = Dict{Vector{Int64}, Function}()  # 创建散射角插值函数字典，键为[入射原子类型, 目标原子类型]
    τFunctions = Dict{Vector{Int64}, Function}()  # 创建飞行时间插值函数字典，键为[入射原子类型, 目标原子类型]
    log_separator()  # 输出日志分隔符
    log_info("Loading θ and τ functions...")  # 记录加载θ和τ函数的信息
    # 遍历所有可能的原子类型对组合
    for type_p in keys(typeDict)  # type_p: 入射原子类型
        for type_t in keys(typeDict)  # type_t: 目标原子类型
            mass_p = typeDict[type_p].mass  # 获取入射原子质量(单位：u)
            mass_t = typeDict[type_t].mass  # 获取目标原子质量(单位：u)
            # 调用θτFunctions函数获取该原子类型对的插值器
            θInterpolation, τInterpolation = θτFunctions(mass_p, mass_t, type_p, type_t, constantsByType, parameters)
            # 创建散射角插值函数：输入能量对数E_p和碰撞参数对数p，输出散射角θ
            θFunctions[[type_p, type_t]] = (E_p, p) -> θInterpolation(E_p, p)
            # 创建飞行时间插值函数：输入能量对数E_p和碰撞参数对数p，输出飞行时间参数τ
            τFunctions[[type_p, type_t]] = (E_p, p) -> τInterpolation(E_p, p)
            # 记录调试信息，显示已加载的原子类型对
            log_debug("  $(parameters.typeDict[type_p].name) → $(parameters.typeDict[type_t].name) loaded")
        end
    end
    log_success("All θ and τ functions initialized")  # 记录成功初始化所有函数
    log_separator()  # 输出日志分隔符
    return θFunctions, τFunctions  # 返回散射角和飞行时间插值函数字典
end

# 计算或加载特定原子类型对的θ和τ插值函数
function θτFunctions(mass_p::Float64, mass_t::Float64, type_p::Int64, type_t::Int64, constantsByType::ConstantsByType, parameters::Parameters)
    E_p_axis = Float64[]  # 初始化能量对数轴数组
    p_axis = Float64[]    # 初始化碰撞参数对数轴数组
    θMatrix = Matrix{Float64}(undef, 0, 0)  # 初始化散射角矩阵
    τMatrix = Matrix{Float64}(undef, 0, 0)  # 初始化飞行时间参数矩阵
    try
        # 尝试从文件加载预计算的θ和τ数据
        E_p_axis, p_axis, θMatrix, τMatrix = LoadθτData(type_p, type_t, parameters)
    catch
        # 如果加载失败，则重新计算θ和τ数据
        EPowerRange = parameters.EPowerRange    # 获取能量对数范围
        pPowerRange = parameters.pPowerRange   # 获取碰撞参数对数范围
        nE = length(EPowerRange)  # 能量轴点数
        np = length(pPowerRange)  # 碰撞参数轴点数
        θMatrix = Array{Float64, 2}(undef, nE, np)  # 分配散射角矩阵内存
        τMatrix = Array{Float64, 2}(undef, nE, np)  # 分配飞行时间参数矩阵内存
        N = length(EPowerRange)  # 总能量点数
        # 使用多线程并行计算θ和τ矩阵，显示进度条
        @showprogress @threads for i in 1:N
            E_p_power = EPowerRange[i]  # 获取当前能量对数
            E_p = 10.0^E_p_power        # 计算实际能量值(单位：eV)
            for (j, p_power) in enumerate(pPowerRange)
                p = 10.0^p_power  # 计算实际碰撞参数值(单位：Å)
                # 调用BCA模块的θτ函数计算散射角和飞行时间参数
                θ, τ = BCA.θτ(E_p, mass_p, mass_t, type_p, type_t, p, constantsByType)
                θMatrix[i, j] = θ  # 存储散射角到矩阵
                τMatrix[i, j] = τ  # 存储飞行时间参数到矩阵
            end
        end
        E_p_axis = collect(EPowerRange)  # 将能量对数范围转换为数组
        p_axis = collect(pPowerRange)    # 将碰撞参数对数范围转换为数组
        # 将计算得到的θ和τ数据保存到文件，供后续使用
        SaveθτData(type_p, type_t, θMatrix, τMatrix, E_p_axis, p_axis, parameters)
    end
    # 创建二维插值函数
    θFunction = interpolate((E_p_axis, p_axis), θMatrix, Gridded(Linear()))  # 散射角线性插值函数
    τFunction = interpolate((E_p_axis, p_axis), τMatrix, Gridded(Linear()))  # 飞行时间参数线性插值函数
    return θFunction, τFunction  # 返回插值函数
end

# =====================================================================
# 新增：模拟器构造函数 - 模块化版本，接受预创建的原子向量
# 功能：通过盒子对象、原子向量和输入网格向量创建模拟器
# 设计理念：将原子创建与模拟器初始化分离，提高代码的模块化和可重用性
# 参数说明：
#   box::Box - 模拟盒子对象，定义模拟空间边界和周期性
#   atoms::Vector{Atom} - 原子向量，包含所有初始原子
#   inputGridVectors::Matrix{Float64} - 输入网格向量，定义网格单元尺寸
#   parameters::Parameters - 模拟参数集合
# =====================================================================
function Simulator(box::Box, atoms::Vector{Atom}, inputGridVectors::Matrix{Float64}, parameters::Parameters)
    # 输入验证
    if size(inputGridVectors) != (3, 3)
        throw(InvalidParameterError("inputGridVectors", size(inputGridVectors), "Must be a 3×3 matrix"))
    end
    det_grid = det(inputGridVectors)
    if det_grid <= 0.0
        throw(InvalidParameterError("inputGridVectors", det_grid, "Determinant must be positive (volume > 0)"))
    end
    if !parameters.is_dynamic_load && isempty(atoms)
        throw(InvalidParameterError("atoms", length(atoms), "Atoms vector cannot be empty in static load mode"))
    end
    
    log_section("Initializing Simulator")  # 记录模拟器初始化开始
    # 调用内部构造函数创建模拟器实例，此时原子和晶格点尚未加载
    simulator = Simulator(box, inputGridVectors, parameters)
    if !parameters.is_dynamic_load  # 如果不是动态加载模式
        LoadAtoms!(simulator, atoms)  # 加载原子和晶格点到模拟器中
    end
    log_success("Simulator initialized.")  # 记录模拟器初始化成功
    return simulator 
end

# =====================================================================
# 新增：原子向量加载函数
# 功能：将预创建的原子向量批量加载到模拟器中，并初始化相关数据结构
# 设计优势：相比逐个原子创建，批量加载提高效率，便于原子向量的重用
# 参数说明：
#   simulator::Simulator - 模拟器对象，原子将被加载到其中
#   atoms::Vector{Atom} - 预创建的原子向量
# =====================================================================
function LoadAtoms!(simulator::Simulator, atoms::Vector{Atom})
    # 遍历原子向量，将每个原子添加到模拟器并创建对应的晶格点
    for atom in atoms
        push!(simulator, atom) 
        latticePoint = LatticePoint(atom)
        push!(simulator, latticePoint)
    end
    # 计算每个网格单元的原子密度（用于非局域能量损失计算）
    for cell in simulator.grid.cells
        cell.atomicDensity = length(cell.atoms) / simulator.grid.cellVolume
    end 
    # 对每个原子施加热扰动（考虑温度效应和非晶化）
    for atom in simulator.atoms
        Pertubation!(atom, simulator)
    end
    # 初始化晶格点环境信息（用于基于环境的DTE计算）
    InitLatticePointEnvronment(simulator)
    log_success("$(simulator.numberOfAtoms) atoms loaded.")  # 记录原子加载成功
end

# =====================================================================
# 新增：根据原胞参数创建原子向量
# 功能：根据原胞基向量、晶格范围和原胞内原子位置生成完整的原子向量
# 设计优势：分离原子创建逻辑，便于原子向量的重用和测试
# 参数说明：
#   parameters::Parameters - 模拟参数，包含原胞定义和晶格信息
# 返回：包含所有晶格原子的向量，按晶格位置顺序排列
# =====================================================================
function CreateAtomsByPrimaryVectors(parameters::Parameters)
    primaryVectors = parameters.primaryVectors  # 原胞基向量
    latticeRanges = parameters.latticeRanges    # 晶格范围
    basis = parameters.basis                    # 原胞内原子基矢坐标
    basisTypes = parameters.basisTypes          # 原胞内原子类型
    # 计算总原子数
    atomNumber = (latticeRanges[1,2] - latticeRanges[1,1]) * (latticeRanges[2,2] - latticeRanges[2,1]) * (latticeRanges[3,2] - latticeRanges[3,1]) * length(basisTypes)
    atoms = Vector{Atom}(undef, atomNumber)  # 预分配原子向量
    n = 1  # 原子向量索引
    
    # 使用进度条显示原子创建过程
    @showprogress desc="Creating atoms ($(atomNumber)): " for x in latticeRanges[1,1]:latticeRanges[1,2]-1
        for y in latticeRanges[2,1]:latticeRanges[2,2]-1    
            for z in latticeRanges[3,1]:latticeRanges[3,2]-1
                for i in 1:length(basisTypes)  # 遍历原胞内所有原子位置
                    reducedCoordinate = Float64[x,y,z] + basis[i, :]  # 计算约化坐标（晶格坐标+基矢偏移）
                    coordinate = primaryVectors' * reducedCoordinate  # 将约化坐标转换为实际空间坐标
                    atoms[n] = Atom(basisTypes[i], coordinate, parameters)  # 创建原子并存入向量
                    n += 1  # 递增索引
                end
            end
        end
    end
    return atoms  # 返回原子向量
end

# =====================================================================
# 修改：模拟器构造函数 - 通过盒子向量和输入网格向量创建模拟器
# 更新：现在使用新的模块化设计，先创建原子向量，再调用新的Simulator构造函数
# 参数说明：
#   boxVectors::Matrix{Float64} - 盒子基向量矩阵，定义模拟空间尺寸
#   inputGridVectors::Matrix{Float64} - 输入网格向量，定义网格单元尺寸
#   parameters::Parameters - 模拟参数集合
# =====================================================================
function Simulator(boxVectors::Matrix{Float64}, inputGridVectors::Matrix{Float64}, parameters::Parameters)
    box = Box(boxVectors)  # 直接通过盒子向量创建模拟盒子
    # 如果不是动态加载模式，则创建原子向量；否则，在动态加载模式下，原子将在需要时加载
    if !parameters.is_dynamic_load
        atoms = CreateAtomsByPrimaryVectors(parameters)
    else
        atoms = Vector{Atom}()  # 动态加载模式下，初始原子向量为空
    end
    # 调用新的Simulator构造函数，传入原子向量
    simulator = Simulator(box, atoms, inputGridVectors, parameters)
    return simulator    
end 

# =====================================================================
# 修改：模拟器构造函数 - 通过盒子尺寸和输入网格向量创建模拟器
# 更新：现在使用新的模块化设计，先创建原子向量，再调用新的Simulator构造函数
# 参数说明：
#   boxSizes::Vector{Int64} - 盒子各维度的原胞数量
#   inputGridVectors::Matrix{Float64} - 输入网格向量，定义网格单元尺寸
#   parameters::Parameters - 模拟参数集合
# =====================================================================
function Simulator(boxSizes::Vector{Int64}, inputGridVectors::Matrix{Float64}, parameters::Parameters)
    # 输入验证
    if length(boxSizes) != 3
        throw(InvalidParameterError("boxSizes", length(boxSizes), "Must be a 3-element vector"))
    end
    if any(s -> s <= 0, boxSizes)
        throw(InvalidParameterError("boxSizes", boxSizes, "All sizes must be positive"))
    end
    
    box = CreateBoxByPrimaryVectors(parameters.primaryVectors, boxSizes)  # 根据原胞矢量和盒子尺寸创建模拟盒子
    # 如果不是动态加载模式，则创建原子向量；否则，在动态加载模式下，原子将在需要时加载
    if !parameters.is_dynamic_load
        atoms = CreateAtomsByPrimaryVectors(parameters)
    else
        atoms = Vector{Atom}()  # 动态加载模式下，初始原子向量为空
    end
    # 调用新的Simulator构造函数，传入原子向量
    simulator = Simulator(box, atoms, inputGridVectors, parameters)
    return simulator    
end

# 参数构造函数：简化版本，使用默认晶体结构
function Parameters(pMax::Float64, vacancyRecoverDistance::Float64, typeDict::Dict{Int64, Element}; kwargs...)
    # 非晶格信息参数
    primaryVectors = [1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0]  # 默认原胞矢量：单位矩阵（简单立方）
    latticeRanges = [0 1; 0 1; 0 1]  # 默认晶格范围：单个原胞
    basis = [0.0 0.0 0.0]  # 默认基矢：单个原子在原点
    basisTypes = [1]  # 默认原子类型：类型1
    # 调用完整参数构造函数，传入默认值和额外关键字参数
    parameters = Parameters(primaryVectors, latticeRanges, basisTypes, basis, pMax, vacancyRecoverDistance, typeDict; kwargs...)
    return parameters  # 返回参数实例
end

# =====================================================================
# 新增：从数据文件加载原子和盒子信息
# 功能：读取数据文件中的原子坐标和类型，创建相应的原子向量和盒子对象
# 设计优势：支持原胞复制，便于从分子动力学模拟快照创建初始结构
# 参数说明：
#   fileName::String - 数据文件名，包含原子坐标和类型信息
#   replicate::Vector{Int64} - 原胞复制因子，默认[1,1,1]表示不复制
# 返回：(box, atoms) 元组，包含创建的盒子对象和原子向量
# =====================================================================
function LoadAtomsAndBoxFromDataFile(fileName::String; replicate::Vector{Int64} = [1,1,1])
    # 从数据文件读取原子信息：盒子边界、原子类型和坐标，并按照replicate参数复制原胞
    xlo, xhi, ylo, yhi, zlo, zhi, types, xs, ys, zs = ReadDate(fileName, replicate)
    # 根据读取的边界创建模拟盒子
    box = Box([xhi-xlo 0.0 0.0; 0.0 yhi-ylo 0.0; 0.0 0.0 zhi-zlo])
    atoms = Vector{Atom}(undef, length(types))  # 预分配原子向量
    # 遍历所有原子数据，创建原子对象
    for (n,(type, x, y, z)) in enumerate(zip(types, xs, ys, zs))
        atoms[n] = Atom(type, [x, y, z], parameters)  # 创建原子并存入向量
    end
    return box, atoms  # 返回盒子对象和原子向量
end

# =====================================================================
# 修改：模拟器构造函数 - 从数据文件加载原子位置
# 更新：使用新的模块化加载函数，支持原胞复制功能
# 参数说明：
#   fileName::String - 数据文件名，包含原子坐标和类型信息
#   inputGridVectors::Matrix{Float64} - 输入网格向量，定义网格单元尺寸
#   parameters::Parameters - 模拟参数集合
#   replicate::Vector{Int64} - 原胞复制因子，默认[1,1,1]表示不复制
# 注意：动态加载模式下不支持从文件加载
# =====================================================================
function Simulator(fileName::String, inputGridVectors::Matrix{Float64}, parameters::Parameters; replicate::Vector{Int64} = [1,1,1])
    if parameters.is_dynamic_load  # 检查是否处于动态加载模式
        throw(SimulationError("Simulator construction", "Simulator from date file is not supported in dynamic load mode."))
    end 
    # 从数据文件加载原子和盒子
    box, atoms = LoadAtomsAndBoxFromDataFile(fileName; replicate=replicate)
    # 调用新的Simulator构造函数，传入加载得到的盒子和原子向量
    simulator = Simulator(box, atoms, inputGridVectors, parameters)
    return simulator
end

# 晶格点构造函数：从原子创建对应的晶格点
function LatticePoint(atom::Atom)
    environment = Vector{Int64}()  # 初始化空的邻近环境向量
    # 创建晶格点实例，复制原子的索引、类型、坐标和单元索引
    return LatticePoint(copy(atom.index), copy(atom.type), 
                        copy(atom.coordinate), atom.cellIndex, environment,
                        atom.index)  # 初始时原子占据对应的晶格点
end

# 确定坐标所在的网格单元
function WhichCell(coordinate::Vector{Float64}, grid::Grid)
    cellIndex = Vector{Int64}(undef, 3)  # 初始化3维单元索引向量
    for d in 1:3  # 遍历三个维度
        # 计算坐标在该维度上的单元索引：坐标/单元尺寸，向下取整
        cellIndex[d] = Int64(floor(coordinate[d] / grid.vectors[d,d])) + 1
        # 边界检查：确保索引在有效范围内
        if cellIndex[d] < 1 
            cellIndex[d] = 1  # 低于下限时设为1
        elseif cellIndex[d] > grid.sizes[d]
            cellIndex[d] = grid.sizes[d]  # 超过上限时设为最大值
        end
    end
    return (cellIndex[1], cellIndex[2], cellIndex[3])  # 返回三元组形式的单元索引
end

# 向模拟器中添加原子的重载push!函数
function push!(simulator::Simulator, atom::Atom)
    atom.index = simulator.maxAtomID + 1  # 分配新的原子索引
    simulator.maxAtomID += 1  # 更新最大原子ID
    push!(simulator.atoms, atom)  # 将原子添加到模拟器的原子列表
    simulator.numberOfAtoms += 1  # 增加原子计数
    cellIndex = WhichCell(atom.coordinate, simulator.grid)  # 确定原子所在的网格单元
    atom.cellIndex = cellIndex  # 设置原子的单元索引
    push!(GetCell(simulator.grid, cellIndex).atoms, atom)  # 将原子添加到对应网格单元的原子列表
end 

# 将晶格点添加到模拟器中，并建立原子与晶格点的双向关联
function push!(simulator::Simulator, latticePoint::LatticePoint)
    push!(simulator.latticePoints, latticePoint)  # 将晶格点添加到模拟器的晶格点列表中
    push!(GetCell(simulator.grid, latticePoint.cellIndex).latticePoints, latticePoint)  # 将晶格点添加到对应网格单元的晶格点列表中
    simulator.atoms[latticePoint.atomIndex].latticePointIndex = latticePoint.index  # 建立原子到晶格点的反向引用
end 

# 将原子添加到网格单元中，并更新原子密度（仅限静态加载模式）
function push!(cell::Cell, atom::Atom, simulator::Simulator)
    atom.cellIndex = cell.index  # 设置原子的网格单元索引
    push!(cell.atoms, atom)  # 将原子添加到网格单元的原子列表中
    if !simulator.parameters.is_dynamic_load  # 仅在静态加载模式下更新原子密度
        cell.atomicDensity = length(cell.atoms) / simulator.grid.cellVolume  # 重新计算网格单元的原子数密度
    end
end 

# 从模拟器中删除原子，处理晶格点关联和状态更新（静态加载版本）
function delete_staticLoad!(simulator::Simulator, atom::Atom)
    originalCell = GetCell(simulator.grid, atom.cellIndex)  # 获取原子当前所在的网格单元
    deleteat!(originalCell.atoms, findfirst(a -> a.index == atom.index, originalCell.atoms))  # 从网格单元原子列表中删除该原子
    simulator.numberOfAtoms -= 1  # 减少模拟器中的原子总数计数
    atom.isAlive = false  # 将原子标记为非活动状态
    if atom.latticePointIndex != -1  # 如果原子占据晶格点位置
        LeaveLatticePoint_staticLoad!(atom, simulator)  # 处理原子离开晶格点的相关操作
    end
end

# 处理原子离开晶格点的过程，创建空位并更新相关数据结构（静态加载版本）
function LeaveLatticePoint_staticLoad!(atom::Atom, simulator::Simulator; isUpdateEnv::Bool = true)
    AddToStore!(atom, simulator)  # 将原子添加到存储中（用于后续恢复）
    latticePoint = simulator.latticePoints[atom.latticePointIndex]  # 获取原子对应的晶格点
    latticePoint.atomIndex = -1  # 将晶格点标记为空位（无原子占据）
    atom.latticePointIndex = -1  # 清除原子的晶格点关联
    vacancy = Atom(latticePoint.type, latticePoint.coordinate, simulator.parameters)  # 创建空位原子对象
    vacancy.index = latticePoint.index  # 设置空位的索引与晶格点一致
    push!(simulator.vacancies, vacancy)  # 将空位添加到模拟器的空位列表中
    push!(GetCell(simulator.grid, latticePoint.cellIndex).vacancies, vacancy)  # 将空位添加到对应网格单元的空位列表中
    if isUpdateEnv && simulator.parameters.isKMC  # 如果需要更新环境且启用了KMC模拟
        DeleteAtomEvents!(simulator, atom)  # 删除原子的KMC事件
        UpdateEvents!(Set(latticePoint.environment), simulator)  # 更新相关晶格点的KMC事件
    end
end

# 从网格单元中删除原子，处理原子密度更新（仅限静态加载模式）
function delete!(cell::Cell, atom::Atom, simulator::Simulator)
    if !atom.isAlive  # 安全检查：确保原子处于活动状态
        error("Atom $(atom.index) is not alive when deleting")  # 抛出错误信息
    end
    #@show atom.index, atom.cellIndex, simulator.nCascade, simulator.nCollisionEvent, atom.coordinate  # 调试输出（注释状态）
    deleteat!(cell.atoms, findfirst(a -> a.index == atom.index, cell.atoms))  # 从网格单元原子列表中删除该原子
    atom.cellIndex = (-1, -1, -1)  # 将原子的网格单元索引设置为无效值
    if !simulator.parameters.is_dynamic_load  # 仅在静态加载模式下更新原子密度
        cell.atomicDensity = length(cell.atoms) / simulator.grid.cellVolume  # 重新计算网格单元的原子数密度
    end
end

# 将原子移动到新位置，处理周期性边界条件和网格单元更新（静态向量版本）
function DisplaceAtom!(atom::Atom, newPosition::SVector{3,Float64}, simulator::Simulator)
    pos = if newPosition isa SVector  # 检查输入是否为静态向量
        [newPosition[1], newPosition[2], newPosition[3]]  # 转换为普通向量以便修改
    else
        copy(newPosition)  # 复制向量以避免修改原始数据
    end
    
    # 处理周期性边界条件和边界反射
    for d in 1:3  # 遍历三个维度
        # 需要适应非周期性边界条件
        if pos[d] < 0  # 如果位置超出下边界
            if simulator.parameters.periodic[d] == false  # 如果该维度是非周期性的
                pos[d] = 0.01  # 反射到边界内部（避免精确边界值）
            else  # 如果该维度是周期性的
                pos[d] += simulator.box.vectors[d,d]  # 应用周期性边界条件（正向穿越）
            end
        elseif pos[d] >= simulator.box.vectors[d,d]  # 如果位置超出上边界
            if simulator.parameters.periodic[d] == false  # 如果该维度是非周期性的
                pos[d] = simulator.box.vectors[d,d] - 0.01  # 反射到边界内部（避免精确边界值）
            else  # 如果该维度是周期性的
                pos[d] -= simulator.box.vectors[d,d]  # 应用周期性边界条件（负向穿越）
            end
        end
    end
    
    SetCoordinate!(atom, pos)  # 设置原子的新坐标
    cellIndex = WhichCell(atom.coordinate, simulator.grid)  # 计算原子新位置所在的网格单元

    if cellIndex != atom.cellIndex  # 如果原子跨越了网格单元边界
        ChangeCell!(atom, cellIndex, simulator)  # 将原子转移到新的网格单元
    end
end

# 将原子移动到新位置的重载函数（普通向量版本）
function DisplaceAtom!(atom::Atom, newPosition::Vector{Float64}, simulator::Simulator)
    DisplaceAtom!(atom, SVector{3,Float64}(newPosition[1], newPosition[2], newPosition[3]), simulator)  # 转换为静态向量并调用主函数
end

# 计算两点间距离的平方，考虑周期性边界条件
function ComputeDistance_squared(coordinate1::Vector{Float64}, coordinate2::Vector{Float64}, crossFlag::NTuple{3, Int8}, box::Box)
    dv = VectorDifference(coordinate1, coordinate2, crossFlag, box)  # 计算考虑周期性边界的向量差
    distance_squared = dv[1]* dv[1] + dv[2]*dv[2] + dv[3]  * dv[3]  # 计算欧几里得距离的平方
    return distance_squared  # 返回距离平方值（避免开方运算提高性能）
end

# 计算两点间的实际距离，考虑周期性边界条件
function ComputeDistance(coordinate1::Vector{Float64}, coordinate2::Vector{Float64}, crossFlag::NTuple{3, Int8}, box::Box)
    return sqrt(ComputeDistance_squared(coordinate1, coordinate2, crossFlag, box))  # 计算实际欧几里得距离
end

# 计算入射原子到目标原子的速度方向距离（投影距离）
function ComputeVDistance(atom_p::Atom, atom_t::Atom, crossFlag::NTuple{3, Int8}, box::Box)
    # v 表示速度方向上的投影
    dv = VectorDifference(atom_p.coordinate, atom_t.coordinate, crossFlag, box)  # 计算考虑周期性的位置向量差
    return dv' * atom_p.velocityDirection  # 返回向量差在入射原子速度方向上的投影（标量积）
end

# 计算两个向量的差，考虑周期性边界穿越
function VectorDifference(v1::Vector{Float64}, v2::Vector{Float64}, crossFlag::NTuple{3, Int8}, box::Box)
    if crossFlag == (Int8(0), Int8(0), Int8(0))  # 如果没有边界穿越
        return v2 - v1  # 直接计算向量差
    end 
    # 考虑周期性边界穿越的向量差计算
    return SVector{3,Float64}(
        v2[1] - v1[1] + crossFlag[1] * box.vectors[1,1],  # x分量：直接差加上边界穿越修正
        v2[2] - v1[2] + crossFlag[2] * box.vectors[2,2],  # y分量：直接差加上边界穿越修正
        v2[3] - v1[3] + crossFlag[3] * box.vectors[3,3]   # z分量：直接差加上边界穿越修正
    )
end

# 计算碰撞参数p和相关几何量，用于二进制碰撞近似计算
function ComputeP!(atom_p::Atom, atom_t::Atom, crossFlag::NTuple{3, Int8}, box::Box)
    dv = VectorDifference(atom_p.coordinate, atom_t.coordinate, crossFlag, box)  # 计算考虑周期性的位置向量差
    t = dot(dv, atom_p.velocityDirection)  # 计算沿速度方向的投影参数t（标量）
    atom_t.pL = t  # 存储投影距离参数（用于非局域能量损失计算）
    
    # 计算碰撞点坐标，考虑周期性边界穿越
    if 1 in crossFlag || -1 in crossFlag  # 如果存在边界穿越
        pPoint_calc = Vector{Float64}(atom_p.coordinate + t * atom_p.velocityDirection)  # 计算初步碰撞点坐标
        for d in 1:3  # 遍历三个维度
            if crossFlag[d] != 0  # 如果该维度有边界穿越
                pPoint_calc[d] -= crossFlag[d] * box.vectors[d,d]  # 修正碰撞点坐标（消除周期性偏移）
            end
        end
    else  # 如果没有边界穿越
        pPoint_calc = atom_p.coordinate + t * atom_p.velocityDirection  # 直接计算碰撞点坐标
    end
    
    atom_t.pPoint = SVector{3,Float64}(pPoint_calc[1], pPoint_calc[2], pPoint_calc[3])  # 存储碰撞点坐标（静态向量）
    pVector_calc = atom_t.pPoint - atom_t.coordinate  # 计算从目标原子指向碰撞点的向量
    atom_t.pVector = SVector{3,Float64}(pVector_calc[1], pVector_calc[2], pVector_calc[3])  # 存储碰撞向量（静态向量）
    p = norm(atom_t.pVector)  # 计算碰撞参数p的大小（目标原子到碰撞路径的垂直距离）
    atom_t.pValue = p  # 存储碰撞参数p的值
    # 需要检查周期性边界条件（注释：可能需要进一步验证）
    return p  # 返回碰撞参数p
end

# 将原子从一个网格单元转移到另一个网格单元
function ChangeCell!(atom::Atom, nextCellIndex::Tuple{Int64, Int64, Int64}, simulator::Simulator)
    originalCell = GetCell(simulator.grid, atom.cellIndex)  # 获取原子当前所在的网格单元
    delete!(originalCell, atom, simulator)  # 从原网格单元中删除原子
    nextCell = GetCell(simulator.grid, nextCellIndex)  # 获取目标网格单元
    push!(nextCell, atom, simulator)  # 将原子添加到目标网格单元
end

# 设置原子的速度方向向量，进行归一化处理
function SetVelocityDirection!(atom::Atom, velocity::SVector{3,Float64})
    n = norm(velocity)  # 计算速度向量的模长
    if isnan(n) || n == Inf || n == 0.0  # 检查无效的模长值（NaN、无穷大或零）
        atom.velocityDirection = SVector{3,Float64}(0.0, 0.0, 0.0)  # 设置速度方向为零向量
    else  # 有效的模长值
        normalized_velocity = velocity / n  # 计算归一化的速度方向向量
        atom.velocityDirection = SVector{3,Float64}(normalized_velocity[1], normalized_velocity[2], normalized_velocity[3])  # 存储归一化后的速度方向
    end
end

# 设置原子速度方向的重载函数（普通向量版本）
function SetVelocityDirection!(atom::Atom, velocity::Vector{Float64})
    SetVelocityDirection!(atom, SVector{3,Float64}(velocity[1], velocity[2], velocity[3]))  # 转换为静态向量并调用主函数
end

# 设置原子的动能，确保能量值为非负
function SetEnergy!(atom::Atom, energy::Float64)
    if energy < 0.0  # 如果能量为负值
        atom.energy = 0.0  # 将能量设置为零（物理约束）
    else  # 能量为非负值
        atom.energy = energy  # 直接设置能量值
    end
end

# 在邻近区域搜索最近的空位，用于空位恢复机制
function GetNeighborVacancy(atom::Atom, simulator::Simulator)
    grid = simulator.grid  # 获取模拟器的网格系统
    cell = GetCell(grid, atom.cellIndex)  # 获取原子当前所在的网格单元
    nearestVacancyDistance_squared = Inf  # 初始化最近空位距离平方为无穷大
    nearestVacancyIndex = -1  # 初始化最近空位索引为无效值
    
    # 遍历所有邻近网格单元
    for neighborCellInfo in cell.neighborCellsInfo
        index, cross = neighborCellInfo.index, neighborCellInfo.cross  # 获取邻近单元索引和边界穿越标志
        neighborCell = GetCell(grid, index)  # 获取邻近网格单元
        
        # 遍历邻近单元中的所有空位
        for vacancy in neighborCell.vacancies
            dr2 = ComputeDistance_squared(atom.coordinate, vacancy.coordinate, cross, simulator.box)  # 计算原子与空位的距离平方
            # 检查距离是否在恢复范围内且是当前找到的最近空位
            if dr2 < simulator.parameters.vacancyRecoverDistance_squared && dr2 < nearestVacancyDistance_squared
                nearestVacancyDistance_squared = dr2  # 更新最近距离平方
                nearestVacancyIndex = vacancy.index  # 更新最近空位索引
            end
        end
    end
    return nearestVacancyIndex  # 返回最近空位的索引（若无合适空位则返回-1）
end

# 停止原子运动并将其恢复到晶格位置（静态加载版本）
function Stop_staticLoad!(atom::Atom, simulator::Simulator)
    SetEnergy!(atom, 0.0)  # 将原子动能设置为0，停止其运动
    SetVelocityDirection!(atom, SVector{3,Float64}([0.0, 0.0, 0.0]))  # 将原子速度方向设置为零向量
    Recover!(atom, simulator)  # 尝试将原子恢复到最近的晶格空位
end

# 将停止的原子恢复到最近的晶格空位
function Recover!(atom::Atom, simulator::Simulator)
    nearestVacancyIndex = GetNeighborVacancy(atom, simulator)  # 查找邻近的空位索引
    if nearestVacancyIndex != -1  # 如果找到邻近空位
        SetOnLatticePoint!(atom, simulator.latticePoints[nearestVacancyIndex], simulator)  # 将原子放置到空位
        # 从模拟器的空位列表中删除该空位
        deleteat!(simulator.vacancies, findfirst(v -> v.index == nearestVacancyIndex, simulator.vacancies))
        # 从所在网格单元的空位列表中删除该空位
        cell = GetCell(simulator.grid, atom.cellIndex)
        deleteat!(cell.vacancies, findfirst(v -> v.index == nearestVacancyIndex, cell.vacancies))
    end
end

# 将原子设置到指定的晶格点上
function SetOnLatticePoint!(atom::Atom, latticePoint::LatticePoint, simulator::Simulator; isUpdateEnv::Bool = true)
    SetEnergy!(atom, 0.0)  # 重置原子能量为0
    SetVelocityDirection!(atom, SVector{3,Float64}([0.0, 0.0, 0.0]))  # 重置速度方向为零
    latticePoint.atomIndex = atom.index  # 将晶格点的原子索引设置为当前原子
    atom.latticePointIndex = latticePoint.index  # 将原子的晶格点索引设置为当前晶格点
    SetCoordinate!(atom, latticePoint.coordinate)  # 将原子坐标设置为晶格点坐标
    # 如果原子存活且不在正确的网格单元中，则移动到正确单元
    if atom.isAlive && atom.cellIndex != latticePoint.cellIndex
        ChangeCell!(atom, latticePoint.cellIndex, simulator)
    elseif !atom.isAlive  # 如果原子不存活（如新恢复的原子）
        atom.isAlive = true  # 标记为存活
        nextCell = GetCell(simulator.grid, latticePoint.cellIndex)  # 获取目标网格单元
        push!(nextCell, atom, simulator)  # 将原子添加到目标单元
    end 
    # 如果启用了KMC模拟且需要更新环境，则更新相关事件
    if simulator.parameters.isKMC && isUpdateEnv
        latticePointIndexs = Set([latticePoint.environment;latticePoint.index])  # 创建包含环境和当前点的索引集合
        UpdateEvents!(latticePointIndexs, simulator)  # 更新这些晶格点相关的事件
    end
    Pertubation!(atom, simulator)  # 对原子位置施加微扰（热振动或非晶化）
end

# 将原子添加到位移存储列表中（用于后续恢复）
function AddToStore!(atom::Atom, simulator::Simulator)
    # 如果模拟器启用了存储功能且原子索引在存储范围内
    if simulator.isStore && atom.index <= simulator.numberOfAtomsWhenStored
        push!(simulator.displacedAtoms, atom.index)  # 将原子索引添加到位移原子列表
    end 
end

# 从位移存储列表中删除原子
function DeleteFromStore!(atom::Atom, simulator::Simulator)
    # 如果模拟器启用了存储功能且原子索引在存储范围内
    if simulator.isStore && atom.index <= simulator.numberOfAtomsWhenStored
        # 从位移原子列表中删除该原子的索引
        deleteat!(simulator.displacedAtoms, findfirst(==(atom.index), simulator.displacedAtoms))
    end
end 

"""
    Restore!(simulator::Simulator)

恢复模拟器到之前保存的状态。

在静态加载模式下，需要先调用 `Save!` 保存状态。此函数会：
- 移除存储后添加的离子
- 清空空位列表
- 将所有位移的原子恢复到原始晶格位置

在动态加载模式下，此函数会清空当前原子和空位列表，重置计数器。

# 参数
- `simulator::Simulator`: 要恢复的模拟器对象

# 抛出异常
- `SimulationError`: 如果网格未初始化（静态模式还需要先保存）

# 示例
```julia
# 静态加载模式
Save!(simulator)
# ... 执行模拟 ...
Restore!(simulator)  # 恢复到保存时的状态

# 动态加载模式
Restore!(simulator)  # 直接清空当前状态
```
"""
function Restore!(simulator::Simulator)
    # 防御性检查
    if isnothing(simulator.grid)
        throw(SimulationError("Restore!", "Simulator grid is not initialized"))
    end
    
    # 静态加载模式需要先保存，动态加载模式不需要
    if !simulator.parameters.is_dynamic_load
        if !simulator.isStore
            throw(SimulationError("Restore!", "Simulator must be saved before restore in static load mode. Call Save!(simulator) first."))
        end
        Restore_staticLoad!(simulator)  # 调用静态加载恢复函数
    else  # 如果是动态加载模式
        # 动态加载模式下，Restore! 只是清空当前状态，不需要先保存
        Restore_dynamicLoad!(simulator)  # 调用动态加载恢复函数
    end
end

# 静态加载模式下的恢复函数
function Restore_staticLoad!(simulator::Simulator)
    # 处理存储后添加的离子（辐照引入的原子）
    for atom in simulator.atoms[simulator.numberOfAtomsWhenStored+1:end]
        # 从它们的网格单元中删除系统中残留的离子
        # 后续将通过设置simulator.atoms = simulator.atoms[1:maxAtomID]来删除simulator.atoms中的离子
        if atom.isAlive  # 如果原子存活
            delete!(GetCell(simulator.grid, atom.cellIndex), atom, simulator)  # 从网格单元中删除
        end
    end
    latticePoints = simulator.latticePoints  # 获取所有晶格点
    # 清空所有网格单元中的空位列表
    for vacancy in simulator.vacancies
        latticePoint = latticePoints[vacancy.index]  # 获取空位对应的晶格点
        cell = GetCell(simulator.grid, latticePoint.cellIndex)  # 获取晶格点所在网格单元
        empty!(cell.vacancies)  # 清空该单元的空位列表
    end
    empty!(simulator.vacancies)  # 清空模拟器的空位列表

    # 将所有位移的原子恢复到它们的原始晶格位置
    for index in Set(simulator.displacedAtoms)
        atom = simulator.atoms[index]  # 获取位移原子
        if atom.index == atom.latticePointIndex  # 如果原子已经在自己的晶格点上
            continue  # 跳过
        end
        latticePoint = simulator.latticePoints[atom.index]  # 获取原子的原始晶格点
        SetOnLatticePoint!(atom, latticePoint, simulator)  # 将原子设置回原始位置
    end
    maxAtomID = simulator.numberOfAtomsWhenStored  # 获取存储时的最大原子ID
    simulator.atoms = simulator.atoms[1:maxAtomID]  # 截断原子列表，移除存储后添加的离子
    simulator.maxAtomID = maxAtomID  # 更新最大原子ID
    simulator.numberOfAtoms  = maxAtomID  # 更新原子数量
    simulator.nCollisionEvent = 0  # 重置碰撞事件计数器
    empty!(simulator.displacedAtoms)  # 清空位移原子列表
end

"""
    Save!(simulator::Simulator)

保存模拟器的当前状态，用于后续恢复。

在静态加载模式下，此函数会检查所有原子是否都在晶格位置上，并记录当前原子数量。
保存后的模拟器可以通过 `Restore!` 恢复到保存时的状态。

# 参数
- `simulator::Simulator`: 要保存的模拟器对象

# 抛出异常
- `SimulationError`: 如果网格未初始化或原子不在晶格位置上

# 示例
```julia
# 保存初始状态
Save!(simulator)

# 执行多次辐照
for i in 1:100
    Irradiation(simulator, 1000.0, ...)
end

# 恢复到初始状态
Restore!(simulator)
```
"""
function Save!(simulator::Simulator)
    # 防御性检查
    if isnothing(simulator.grid)
        throw(SimulationError("Save!", "Simulator grid is not initialized"))
    end
    if isempty(simulator.atoms)
        throw(SimulationError("Save!", "Cannot save simulator with no atoms"))
    end
    
    # 检查所有原子是否都在晶格位置上
    invalid_atoms = Vector{Int64}()
    for atom in simulator.atoms
        if atom.latticePointIndex == -1  # 如果原子不在晶格上
            push!(invalid_atoms, atom.index)
        end
    end
    
    if !isempty(invalid_atoms)
        throw(SimulationError("Save!", "Atoms not on lattice when stored: $(invalid_atoms[1:min(10, length(invalid_atoms))])" * 
            (length(invalid_atoms) > 10 ? " ... (total: $(length(invalid_atoms)))" : "")))
    end
    
    simulator.isStore = true  # 设置存储标志为true
    simulator.numberOfAtomsWhenStored = simulator.numberOfAtoms  # 记录存储时的原子数量
end 

# 获取指定晶格点的环境晶格点（邻近晶格点）
function GetEnvironmentLatticePoints(latticePoint::LatticePoint, simulator::Simulator)
    cellIndex = latticePoint.cellIndex  # 获取晶格点所在网格单元索引
    theCell = GetCell(simulator.grid, cellIndex)  # 获取晶格点所在网格单元
    grid = simulator.grid  # 获取模拟器网格
    cut_squared = simulator.environmentCut^2  # 获取环境截断距离的平方（用于距离比较优化）
    box = simulator.box  # 获取模拟盒子
    environmentLatticePointsIndex = Vector{Int64}()  # 初始化环境晶格点索引列表
    dVectors = Vector{Vector{Float64}}()  # 初始化距离向量列表（用于排序）
    # 遍历所有邻近网格单元
    for neighborCellInfo in theCell.neighborCellsInfo
        index, cross = neighborCellInfo.index, neighborCellInfo.cross  # 获取邻近单元索引和边界穿越标志
        cell = GetCell(grid, index)  # 获取邻近网格单元
        latticePoints = cell.latticePoints  # 获取邻近单元的晶格点
        # 遍历邻近单元中的所有晶格点
        for neighborLatticePoint in latticePoints
            neighborLatticePointIndex = neighborLatticePoint.index  # 获取邻近晶格点索引
            # 检查距离是否在截断范围内且不是自身
            if ComputeDistance_squared(latticePoint.coordinate, neighborLatticePoint.coordinate, cross, box) <= cut_squared && neighborLatticePointIndex != latticePoint.index
                # 计算距离向量并存储
                push!(dVectors, VectorDifference(latticePoint.coordinate, neighborLatticePoint.coordinate, cross, box))
                push!(environmentLatticePointsIndex, neighborLatticePointIndex)  # 存储邻近晶格点索引
            end
        end
    end
    # 按距离向量的x、y、z分量对索引进行排序（确保环境描述的一致性）
    sorted_indices = sortperm(dVectors, by = v -> (v[1], v[2], v[3]))
    environmentLatticePointsIndex = environmentLatticePointsIndex[sorted_indices]
    
    return environmentLatticePointsIndex  # 返回排序后的环境晶格点索引列表
end

# 初始化所有晶格点的环境信息
function InitLatticePointEnvronment(simulator::Simulator)
    # 如果DTE模式不是1（直接）或4（自定义），则需要初始化环境
    if simulator.parameters.DTEMode != 1 && simulator.parameters.DTEMode != 4
        log_info("🌐 Initializing lattice point environment...\n")  # 输出初始化信息
        # 为每个晶格点计算环境
        for latticePoint in simulator.latticePoints
            latticePoint.environment = GetEnvironmentLatticePoints(latticePoint, simulator)
        end
    end
    # 注释：simulator.environmentLength应该从DTEDict中获取（未来实现）
end

# 计算晶格点的环境索引（用于DTE查找）
function GetEnvironmentIndex(latticePoint::LatticePoint, simulator::Simulator)
    environment = latticePoint.environment  # 获取晶格点的环境列表
    latticePoints = simulator.latticePoints  # 获取所有晶格点
    index = 0  # 初始化环境索引
    
    # 使用二进制编码表示环境占据状态
    for i in 1:length(environment)
        if latticePoints[environment[i]].atomIndex != -1  # 如果环境位置被原子占据
            index += 2^(i-1)  # 设置对应的二进制位
        end
    end
    return index + 1  # 返回索引（加1是因为Julia数组索引从1开始）
end 

# 生成高斯分布的随机位移（用于热振动）
function GaussianDeltaX(sigma::Float64, simulator::Union{Simulator, Nothing}=nothing)
    rng = if isnothing(simulator)
        Main.THREAD_RNG[Threads.threadid()]  # 向后兼容：如果没有 simulator，使用全局的
    else
        get_thread_rng(simulator)  # 使用 simulator 中的 RNG
    end
    return randn(rng) * sigma  # 使用线程安全的RNG生成高斯随机数
end

# 对原子位置施加微扰（热振动或非晶化）
function Pertubation!(atom::Atom, simulator::Simulator)
    if simulator.parameters.isAmorphous  # 如果是非晶材料
        rng = get_thread_rng(simulator)  # 获取线程安全的随机数生成器
        # 在网格单元内随机分布原子位置
        atom.coordinate .= GetCell(simulator.grid, atom.cellIndex).ranges[:,1] .+ [rand(rng) * simulator.grid.vectors[d, d] for d in 1:3]
    else  # 如果是晶体材料
        ah = simulator.parameters.amorphousHeight  # 获取非晶区域高度阈值
        if atom.coordinate[3] > ah  # 如果原子在非晶区域以上
            rng = get_thread_rng(simulator)  # 获取线程安全的随机数生成器
            cell = GetCell(simulator.grid, atom.cellIndex)  # 获取原子所在网格单元
            # 在x和y方向随机分布，z方向在非晶区域内随机分布
            atom.coordinate[1] = cell.ranges[1,1] + rand(rng) * simulator.grid.vectors[1, 1]
            atom.coordinate[2] = cell.ranges[2,1] + rand(rng) * simulator.grid.vectors[2, 2]
            base = cell.ranges[3,1] > ah ? cell.ranges[3,1] : ah  # 计算z方向的基准位置
            latticeTop = simulator.parameters.primaryVectors[3,3] * simulator.parameters.latticeRanges[3,2]  # 计算晶格顶部位置
            top = cell.ranges[3,2] < latticeTop ? cell.ranges[3,2] : latticeTop  # 计算z方向的上限位置
            atom.coordinate[3] = base + rand(rng) * (top - base)  # 在z方向随机分布
        else  # 如果原子在晶体区域
            if simulator.parameters.temperature > 0.0  # 如果温度大于0
                # 对所有三个方向施加热振动
                for d in 1:3
                    atom.coordinate[d] += GaussianDeltaX(simulator.constantsByType.sigma[atom.type], simulator)
                end
            end
        end
    end
end

# 设置原子坐标
function SetCoordinate!(atom::Atom, coordinate::Vector{Float64})
    atom.coordinate .= coordinate  # 使用广播赋值更新原子坐标
end 

# 根据温度计算热振动的均方根位移（德拜模型）
function TemperatureToSigma(T::Float64, θ_D::Float64, m_rel::Float64; atol=1e-10, rtol=1e-8)
    if T == 0.0  # 如果温度为0K
        log_debug("Temperature is 0 K")  # 输出调试信息
        return 0  # 返回零位移
    end
    # 定义物理常数
    ħ   = 1.054_571_817e-34      # 约化普朗克常数，单位：J·s
    kB  = 1.380_649_000e-23      # 玻尔兹曼常数，单位：J/K
    amu = 1.660_539_066_60e-27   # 原子质量单位，单位：kg

    M = m_rel * amu  # 计算实际质量，单位：kg
    y_max = θ_D / T  # 计算积分上限

    # 对x/(e^x-1)进行积分（德拜模型的关键积分）
    integrand(x) = x / (exp(x) - 1)  # 定义被积函数
    I, _ = quadgk(integrand, 0.0, y_max; atol, rtol)  # 数值积分

    # 计算均方位移的方差
    σ2 = 3 * ħ^2 / (M * kB * θ_D) * (0.25 + (T/θ_D)^2 * I)
    σ  = sqrt(σ2) * 1e10         # 从米转换为埃（Å）

    return σ  # 返回均方根位移
end

# =====================================================================
# 统一接口函数
# =====================================================================

"""
    LeaveLatticePoint!(atom::Atom, simulator::Simulator; isUpdateEnv::Bool=true)

处理原子离开晶格点的过程的统一接口。

根据 `simulator.parameters.is_dynamic_load` 自动选择静态或动态加载实现。

# 参数
- `atom::Atom`: 要离开晶格点的原子
- `simulator::Simulator`: 模拟器对象
- `isUpdateEnv::Bool=true`: 是否更新环境（用于KMC模拟）
"""
function LeaveLatticePoint!(atom::Atom, simulator::Simulator; isUpdateEnv::Bool=true)
    if simulator.parameters.is_dynamic_load
        LeaveLatticePoint_dynamicLoad!(atom, simulator; isUpdateEnv=isUpdateEnv)
    else
        LeaveLatticePoint_staticLoad!(atom, simulator; isUpdateEnv=isUpdateEnv)
    end
end

"""
    Stop!(atom::Atom, simulator::Simulator)

停止原子运动并将其恢复到晶格位置的统一接口。

根据 `simulator.parameters.is_dynamic_load` 自动选择静态或动态加载实现。

# 参数
- `atom::Atom`: 要停止的原子
- `simulator::Simulator`: 模拟器对象
"""
function Stop!(atom::Atom, simulator::Simulator)
    if simulator.parameters.is_dynamic_load
        Stop_dynamicLoad!(atom, simulator)
    else
        Stop_staticLoad!(atom, simulator)
    end
end

"""
    delete!(simulator::Simulator, atom::Atom)

从模拟器中删除原子的统一接口。

根据 `simulator.parameters.is_dynamic_load` 自动选择静态或动态加载实现。

# 参数
- `simulator::Simulator`: 模拟器对象
- `atom::Atom`: 要删除的原子
"""
function Base.delete!(simulator::Simulator, atom::Atom)
    if simulator.parameters.is_dynamic_load
        delete_dynamicLoad!(simulator, atom)
    else
        delete_staticLoad!(simulator, atom)
    end
end