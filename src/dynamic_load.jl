# 导入StaticArrays包用于高效的静态数组操作
using StaticArrays

# 正交晶格原子计算函数 - 针对正交晶系的优化版本
function ComputeLatticeAtoms_Orthogonal!(cell::Cell, simulator::Simulator)
    # 获取模拟参数和晶格信息
    parameters = simulator.parameters
    primaryVectors = parameters.primaryVectors  # 原胞基向量矩阵(3×3)
    latticeRanges = parameters.latticeRanges    # 晶格索引范围矩阵(3×2)
    basisTypes = parameters.basisTypes          # 原胞内原子类型向量
    basis = parameters.basis                    # 原胞内原子分数坐标矩阵

    # 如果单元晶格范围未保存，则计算并保存晶格索引范围
    if !cell.isSavedLatticeRange
        # 针对正交情况的优化计算
        # 注释：a1, a2, a3 = primaryVectors[1,1], primaryVectors[2,2], primaryVectors[3,3]

        # 遍历三个维度计算晶格索引范围
        for d in 1:3
            # 计算该维度最小晶格索引：物理坐标转换为晶格坐标，并与全局范围取最大值
            cell.latticeRanges[d, 1] = max(floor(Int, cell.ranges[d, 1] / primaryVectors[d, d]), latticeRanges[d, 1])
            # 计算该维度最大晶格索引：物理坐标转换为晶格坐标，并与全局范围取最小值（减1确保在范围内）
            cell.latticeRanges[d, 2] = min(floor(Int, cell.ranges[d, 2] / primaryVectors[d, d]), latticeRanges[d, 2] - 1)
        end
        cell.isSavedLatticeRange = true  # 标记晶格范围已保存
    end

    # 使用线程安全的工作缓冲区坐标，避免内存分配
    # 注释：coordinate = Vector{Float64}(undef, 3)  # 原始分配方式
    coordinate = simulator.workBuffers.coordinates[Threads.threadid()]  # 线程本地坐标缓冲区
    indexInCell = 0  # 初始化单元内原子索引计数器

    # 三重循环遍历晶格索引范围生成原子
    for x in cell.latticeRanges[1, 1]:cell.latticeRanges[1, 2]
        for y in cell.latticeRanges[2, 1]:cell.latticeRanges[2, 2]
            for z in cell.latticeRanges[3, 1]:cell.latticeRanges[3, 2]
                # 遍历原胞内所有原子位置
                for i in 1:length(basisTypes)
                    indexInCell += 1  # 递增单元内原子索引

                    # 对于非晶材料，跳过被空位占据的位置
                    if !simulator.parameters.isAmorphous && any(v -> v.indexInCell == indexInCell, cell.vacancies)
                        continue  # 该位置被空位占据，跳过原子生成
                    end

                    # 计算原子的物理坐标：晶格索引 + 分数坐标，乘以基向量
                    coordinate[1] = primaryVectors[1, 1] * (x + basis[i, 1])  # x坐标 = a₁*(x + basis_x)
                    coordinate[2] = primaryVectors[2, 2] * (y + basis[i, 2])  # y坐标 = a₂*(y + basis_y)
                    coordinate[3] = primaryVectors[3, 3] * (z + basis[i, 3])  # z坐标 = a₃*(z + basis_z)

                    # 边界检查：仅对边界晶格点进行精确位置验证
                    if x == cell.latticeRanges[1, 1] || x == cell.latticeRanges[1, 2] ||
                       y == cell.latticeRanges[2, 1] || y == cell.latticeRanges[2, 2] ||
                       z == cell.latticeRanges[3, 1] || z == cell.latticeRanges[3, 2]
                        # 检查原子是否确实在当前单元边界内
                        if (coordinate[1] < cell.ranges[1, 1] || coordinate[1] >= cell.ranges[1, 2] ||
                            coordinate[2] < cell.ranges[2, 1] || coordinate[2] >= cell.ranges[2, 2] ||
                            coordinate[3] < cell.ranges[3, 1] || coordinate[3] >= cell.ranges[3, 2])
                            continue  # 原子超出单元边界，跳过
                        end
                    end

                    # 创建原子对象
                    atom = Atom(basisTypes[i], copy(coordinate), parameters)  # 复制坐标创建新原子
                    # 保存理想晶格位置（用于空位恢复）
                    atom.latticeCoordinate = SVector{3,Float64}(atom.coordinate[1], atom.coordinate[2], atom.coordinate[3])
                    atom.cellIndex = cell.index      # 设置原子所在单元索引
                    atom.index = 0                   # 临时索引，将在后续分配全局索引
                    atom.isNewlyLoaded = true        # 标记为新加载的晶格原子
                    atom.indexInCell = indexInCell   # 设置原子在单元内的局部索引

                    Pertubation!(atom, simulator)    # 应用热扰动或非晶化位移

                    push!(cell.latticeAtoms, atom)   # 将原子添加到单元的晶格原子列表

                    # 调试模式下的额外处理（已注释）
                    #if simulator.parameters.debugMode == true
                    #    push!(simulator.debugAtoms, atom)
                    #end
                end
            end
        end
    end

    # 对于非晶材料，调整晶格原子数量以排除空位
    if simulator.parameters.isAmorphous
        resize!(cell.latticeAtoms, length(cell.latticeAtoms) - length(cell.vacancies))
    end
end

# 通用晶格原子计算函数 - 针对非正交晶系（当前实现有缺陷）
function ComputeLatticeAtoms_General!(cell::Cell, simulator::Simulator)  # 注释：this is wrong!!
    parameters = simulator.parameters
    primaryVectors = parameters.primaryVectors      # 原胞基向量
    latticeRanges = parameters.latticeRanges        # 全局晶格索引范围
    basisTypes = parameters.basisTypes              # 原胞原子类型
    basis = parameters.basis                        # 原胞原子分数坐标

    # 如果单元晶格范围未保存，则计算晶格索引范围
    if !cell.isSavedLatticeRange
        # 非正交情况的原始计算方法
        primaryVectors_INV = parameters.primaryVectors_INV  # 原胞基向量的逆矩阵

        # 创建单元8个顶点的矩阵表示
        vertexMatrix = Matrix{Float64}(undef, 8, 3)
        vertexMatrix[1, :] = [ranges[1, 1], ranges[2, 1], ranges[3, 1]]  # 顶点1: (min_x, min_y, min_z)
        vertexMatrix[2, :] = [ranges[1, 1], ranges[2, 1], ranges[3, 2]]  # 顶点2: (min_x, min_y, max_z)
        vertexMatrix[3, :] = [ranges[1, 2], ranges[2, 1], ranges[3, 1]]  # 顶点3: (max_x, min_y, min_z)
        vertexMatrix[4, :] = [ranges[1, 2], ranges[2, 1], ranges[3, 2]]  # 顶点4: (max_x, min_y, max_z)
        vertexMatrix[5, :] = [ranges[1, 1], ranges[2, 2], ranges[3, 1]]  # 顶点5: (min_x, max_y, min_z)
        vertexMatrix[6, :] = [ranges[1, 1], ranges[2, 2], ranges[3, 2]]  # 顶点6: (min_x, max_y, max_z)
        vertexMatrix[7, :] = [ranges[1, 2], ranges[2, 2], ranges[3, 1]]  # 顶点7: (max_x, max_y, min_z)
        vertexMatrix[8, :] = [ranges[1, 2], ranges[2, 2], ranges[3, 2]]  # 顶点8: (max_x, max_y, max_z)

        # 将顶点坐标转换为分数坐标
        nfrac = vertexMatrix * primaryVectors_INV

        # 计算需要覆盖的晶格索引范围（扩大范围确保完全覆盖）
        nmin = [floor(Int, minimum(nfrac[:, 1]) - 1),  # x方向最小索引，减1确保覆盖
            floor(Int, minimum(nfrac[:, 2]) - 1),  # y方向最小索引
            floor(Int, minimum(nfrac[:, 3]) - 1)]  # z方向最小索引
        nmax = [ceil(Int, maximum(nfrac[:, 1]) + 1),   # x方向最大索引，加1确保覆盖
            ceil(Int, maximum(nfrac[:, 2]) + 1),   # y方向最大索引
            ceil(Int, maximum(nfrac[:, 3]) + 1)]   # z方向最大索引

        # 将计算的范围与全局范围进行裁剪
        for d in 1:3
            cell.latticeRanges[d, 1] = max(nmin[d], latticeRanges[d, 1])  # 取最大值确保在全局范围内
            cell.latticeRanges[d, 2] = min(nmax[d], latticeRanges[d, 2])  # 取最小值确保在全局范围内
        end
        cell.isSavedLatticeRange = true  # 标记晶格范围已保存
    else
        # 如果范围已保存，直接使用保存的值
        n1, n2, n3 = cell.latticeRanges
    end

    # 获取单元边界用于坐标过滤
    ranges = cell.ranges
    x_min, x_max = ranges[1, 1], ranges[1, 2]  # x方向物理边界
    y_min, y_max = ranges[2, 1], ranges[2, 2]  # y方向物理边界
    z_min, z_max = ranges[3, 1], ranges[3, 2]  # z方向物理边界

    # 非正交情况的原子生成循环
    for x in n1, y in n2, z in n3, i in 1:length(basisTypes)
        # 计算原子的物理坐标：晶格向量乘以（晶格索引 + 分数坐标）
        coordinate = primaryVectors' * (Float64[x, y, z] + basis[i, :])

        # 检查原子是否在单元边界内
        if (coordinate[1] >= x_min && coordinate[1] <= x_max &&
            coordinate[2] >= y_min && coordinate[2] <= y_max &&
            coordinate[3] >= z_min && coordinate[3] <= z_max)

            # 检查该位置是否被空位占据
            for vacancy in cell.vacancies
                if ComputeDistance_squared(coordinate, vacancy.coordinate, (Int8(0), Int8(0), Int8(0)), simulator.box) < 1E-10
                    continue  # 位置被空位占据，跳过原子生成
                end
            end

            # 创建原子对象
            atom = Atom(basisTypes[i], coordinate, parameters)
            atom.latticeCoordinate = SVector{3,Float64}(atom.coordinate[1], atom.coordinate[2], atom.coordinate[3])
            atom.cellIndex = cell.index        # 设置原子所在单元索引
            atom.index = 0                     # 临时值，将在后续分配
            atom.isNewlyLoaded = true          # 标记为新加载的晶格原子
            Pertubation!(atom, simulator)      # 应用热扰动
            push!(cell.latticeAtoms, atom)     # 将原子添加到晶格原子列表
        end
    end
end

# 加载单元原子函数 - 主入口点
function LoadCellAtoms!(cell::Cell, simulator::Simulator)
    # 如果单元未加载，则加载晶格原子
    if !cell.isLoaded
        # 根据原胞是否正交选择相应的计算函数
        if simulator.parameters.isPrimaryVectorOrthogonal
            ComputeLatticeAtoms_Orthogonal!(cell, simulator)  # 使用正交优化版本
        else
            ComputeLatticeAtoms_General!(cell, simulator)     # 使用通用版本（当前工作不正确）
        end
        # 计算并更新单元的原子密度
        cell.atomicDensity = length(cell.latticeAtoms) / simulator.grid.cellVolume
        cell.isLoaded = true  # 标记单元为已加载状态
    end
end

# 动态加载模式下从邻近单元获取碰撞目标函数
function GetTargetsFromNeighbor_dynamicLoad(atom::Atom, cell::Cell, filterIndexes::Vector{Int64}, simulator::Simulator)
    grid = simulator.grid      # 模拟网格系统
    box = simulator.box        # 模拟盒子
    targets = Vector{Atom}()   # 初始化目标原子向量
    pMax = simulator.parameters.pMax  # 最大碰撞参数
    nthreads = Threads.nthreads()     # 获取可用线程数

    neighborCellsInfo = cell.neighborCellsInfo  # 邻近单元信息数组(3×3×3)

    # 初始化标志数组用于跟踪邻近单元加载状态
    AlreadyLoadedFlags = [true for _ in 1:27]    # 27个邻近单元的加载状态标志
    infiniteFlag_tls = [true for _ in 1:27]      # 线程本地无限标志数组

    # 获取线程候选目标缓冲区并清空
    threadCandidates = simulator.workBuffers.threadCandidates
    for tc in threadCandidates
        empty!(tc)  # 清空每个线程的候选缓冲区
    end

    # 预加载邻近单元以避免竞态条件（多线程安全）
    if nthreads >= 1
        for n in 1:length(neighborCellsInfo)
            neighborCellInfo = neighborCellsInfo[n]
            index = neighborCellInfo.index
            GetCell(grid, index)  # 预加载单元到网格字典中
        end
    end

    # 并行处理所有邻近单元
    @threads :static for n in 1:length(neighborCellsInfo)
        neighborCellInfo = neighborCellsInfo[n]  # 获取当前邻近单元信息
        cross = neighborCellInfo.cross           # 边界穿越标志

        # 检查非周期性边界条件
        nonPeriodicFlag = false
        for d in 1:3
            if cross[d] != 0 && !simulator.parameters.periodic[d]
                nonPeriodicFlag = true  # 在非周期性维度上穿越边界
                break
            end
        end
        if nonPeriodicFlag
            continue  # 跳过非周期性边界的穿越情况
        end

        # 获取当前线程的候选缓冲区
        buf = threadCandidates[Threads.threadid()]
        index = neighborCellInfo.index
        neighborCell = GetCell(grid, index)  # 获取邻近单元

        # 如果邻近单元已被探索，跳过处理
        if neighborCell.isExplored
            continue
        end

        # 更新加载状态标志
        if !neighborCell.isLoaded
            AlreadyLoadedFlags[n] = false  # 标记该单元需要加载
        end

        # 加载邻近单元的原子
        LoadCellAtoms!(neighborCell, simulator)
        neighborCell.isExplored = true    # 标记单元为已探索
        infiniteFlag_tls[n] = false       # 清除无限标志（找到有效单元）

        # 处理邻近单元中的所有原子（常规原子 + 晶格原子）
        na = length(neighborCell.atoms)
        for n in 1:na+length(neighborCell.latticeAtoms)
            # 注释：for neighborAtom in [neighborCell.atoms; neighborCell.latticeAtoms]

            # 选择当前原子：前na个为常规原子，之后为晶格原子
            if n <= na
                neighborAtom = neighborCell.atoms[n]
            else
                neighborAtom = neighborCell.latticeAtoms[n-na]
            end

            # 过滤条件：跳过自身和已过滤的原子
            if neighborAtom.index == atom.index || neighborAtom.index in filterIndexes
                continue
            end

            # 计算速度方向距离，确保原子在运动方向上
            if ComputeVDistance(atom, neighborAtom, neighborCellInfo.cross, box) > 0
                # 计算碰撞参数p
                p = ComputeP!(atom, neighborAtom, neighborCellInfo.cross, box)

                # 如果碰撞参数超过最大值，跳过该原子
                if p >= pMax
                    continue
                end

                # 将候选原子添加到线程缓冲区
                push!(buf, neighborAtom)
            end
        end
    end

    # 合并所有线程的候选目标
    candidateTargets = simulator.workBuffers.candidateTargets
    empty!(candidateTargets)  # 清空全局候选列表
    for tc in threadCandidates
        append!(candidateTargets, tc)  # 合并所有线程的候选
    end

    # 计算整体无限标志：如果所有邻近单元都导致无限飞行
    infiniteFlag = reduce(&, infiniteFlag_tls)

    # 处理已探索的单元和加载状态
    for (neighborCellInfo, flag) in zip(neighborCellsInfo, AlreadyLoadedFlags)
        idx = neighborCellInfo.index
        cell = GetCell(simulator.grid, idx)
        push!(simulator.exploredCells, cell)  # 添加到已探索单元列表

        # 如果单元是新加载的，为晶格原子分配临时负索引
        if flag == false
            for atom in cell.latticeAtoms
                simulator.minLatticeAtomID -= 1  # 递减最小晶格原子ID
                atom.index = simulator.minLatticeAtomID  # 分配临时负索引
            end
        end
    end

    # 如果没有候选目标，返回空结果
    if isempty(candidateTargets)
        return (targets, infiniteFlag)
    end

    # 找到路径长度最小的最近目标
    _, minIdx = findmin(t -> t.pL, candidateTargets)
    nearestTarget = candidateTargets[minIdx]
    push!(targets, nearestTarget)  # 添加最近目标到结果列表

    # 检查其他候选目标是否满足同时碰撞条件
    for candidateTarget in candidateTargets
        if candidateTarget.index == nearestTarget.index
            continue  # 跳过最近目标自身
        end
        # 应用同时碰撞准则检查
        if SimultaneousCriteria(candidateTarget, nearestTarget, simulator)
            push!(targets, candidateTarget)  # 添加满足条件的候选目标
        end
    end

    return (targets, infiniteFlag)  # 返回目标列表和无限标志
end

# 判断候选目标原子是否满足同时碰撞条件的函数
# 根据碰撞动力学理论，当多个目标原子与入射原子的碰撞参数接近时，可能发生同时碰撞
function SimultaneousCriteria(candidateTarget::Atom, nearestTarget::Atom, simulator::Simulator)
    # 计算候选目标与最近目标在飞行路径上的距离差 (单位：Å)
    # 这个差值反映了两个目标原子相对于入射原子轨迹的位置差异
    deltaPL = candidateTarget.pL - nearestTarget.pL

    # 准则1：检查路径长度差是否超过类型相关的最大允许值
    # qMax 是基于原子半径的碰撞参数限制，确保碰撞在物理合理范围内
    if deltaPL > simulator.constantsByType.qMax[[candidateTarget.type, nearestTarget.type]]
        return false  # 距离差过大，不符合同时碰撞条件

    # 准则2：检查最近目标的碰撞参数平方与距离差平方之和是否超过pMax²
    # 这确保了两个碰撞事件在几何上不会超出最大碰撞距离
    elseif nearestTarget.pValue * nearestTarget.pValue + deltaPL * deltaPL > simulator.parameters.pMax_squared
        return false  # 几何约束不满足

    # 准则3：检查候选目标的碰撞参数平方与距离差平方之和是否超过pMax²
    # 从候选目标的角度验证几何约束
    elseif candidateTarget.pValue * candidateTarget.pValue + deltaPL * deltaPL > simulator.parameters.pMax_squared
        return false  # 几何约束不满足
    end

    # 所有准则都满足，返回true表示可以发生同时碰撞
    return true
end

# 动态加载模式下的碰撞处理函数
# 处理入射原子与多个目标原子的同时碰撞，包括能量动量守恒和原子位移
function Collision_dynamicLoad!(atom_p::Atom, atoms_t::Vector{Atom}, simulator::Simulator)
    # 获取目标原子数量
    N_t = length(atoms_t)
    grid = simulator.grid

    # 获取碰撞参数计算缓冲区，避免重复内存分配
    buffers = simulator.workBuffers.collisionParames
    # 确保缓冲区容量足够存储当前碰撞的所有参数
    EnsureCollisionCapacity!(buffers, N_t)

    # 创建缓冲区的视图，用于高效访问碰撞参数数组
    tanφList = @view buffers.tanφList[1:N_t]  # 入射原子散射角正切值数组
    tanψList = @view buffers.tanψList[1:N_t]  # 目标原子散射角正切值数组
    E_tList = @view buffers.E_tList[1:N_t]    # 传递给各目标原子的能量数组(单位：eV)
    x_pList = @view buffers.x_pList[1:N_t]    # 入射原子位移分量数组(单位：Å)
    x_tList = @view buffers.x_tList[1:N_t]    # 目标原子位移分量数组(单位：Å)
    Q_locList = @view buffers.Q_locList[1:N_t] # 局域能量损失数组(单位：eV)

    # 使用第一个目标原子作为参考计算平均参数
    atom_t = atoms_t[1]
    pL = atom_t.pL      # 碰撞路径长度(单位：Å)
    pPoint = atom_t.pPoint  # 碰撞点坐标

    # 修正路径长度，减去入射原子已经飞行的空路径
    pL -= atom_p.emptyPath

    # 获取系统的均匀原子数密度(单位：原子/Å³)
    N = simulator.uniformDensity

    # 调试代码段（已注释）
    # if simulator.nCollisionEvent < 3
    #     @show atom_p.emptyPath
    #     @show pL
    #     @show N
    #     @show simulator.uniformDensity
    # elseif simulator.nCollisionEvent == 3
    #     exit()
    # end

    # 计算非局域电子阻止能量损失（基于Lindhard理论）
    # Q_nl_v 表示沿碰撞路径的连续电子激发能量损失
    Q_nl_v = Q_nl(atom_p.energy, atom_p.mass, atom_t.mass, atom_p.type, atom_t.type,
        pL, N, simulator.constantsByType)

    # 从入射原子能量中减去非局域能量损失
    atom_p.energy -= Q_nl_v

    # 调试代码（已注释）：跟踪特定类型原子的能量损失
    # if atom_p.type == 2
    #     global Q_loss += Q_nl_v  # debug 
    # end

    # 能量修正：如果能量接近停止能量但仍有动能，设置为略高于停止能量
    # 避免数值不稳定和零能量情况
    if atom_p.energy < 0.1 && atom_p.energy + Q_nl_v >= 0.1
        atom_p.energy = 0.11  # 设置为略高于停止能量
    end

    # 初始化总动量向量，用于动量守恒计算
    momentum = @SVector [0.0, 0.0, 0.0]

    # 遍历所有目标原子，计算碰撞参数
    for (i, atom_t) in enumerate(atoms_t)
        p = atom_t.pValue  # 当前目标原子的碰撞参数(单位：Å)

        # 计算二元碰撞近似下的所有碰撞参数
        # 包括散射角、能量转移、位移距离等
        tanφList[i], tanψList[i], E_tList[i], x_pList[i], x_tList[i], Q_locList[i] = CollisionParams(
            atom_p.energy, atom_p.mass, atom_t.mass, atom_p.type, atom_t.type, p, simulator.constantsByType,
            simulator.θFunctions[[atom_p.type, atom_t.type]], simulator.τFunctions[[atom_p.type, atom_t.type]])

        # 计算目标原子碰撞后的速度方向
        # 基于经典碰撞动力学：v_t = - (p_vector/|p|) * tanψ + v_p
        if atom_t.pValue != 0
            velocityDirectionTmp = -atom_t.pVector / atom_t.pValue * tanψList[i] + atom_p.velocityDirection
        else
            # 零碰撞参数情况（对心碰撞），保持入射原子方向
            velocityDirectionTmp = atom_p.velocityDirection
        end

        # 设置目标原子的速度方向
        SetVelocityDirection!(atom_t, velocityDirectionTmp)

        # 累加目标原子的动量贡献
        # 动量 = 质量 × 速度 = sqrt(2 * mass * energy) × direction
        momentum += sqrt(2 * atom_t.mass * E_tList[i]) * atom_t.velocityDirection
    end

    # 计算入射原子碰撞后的动量（动量守恒）
    # p_momentum = 初始动量 - 传递给目标原子的总动量
    pMomentum = sqrt(2 * atom_p.mass * atom_p.energy) * atom_p.velocityDirection - momentum

    # 计算入射原子碰撞后的速度
    pVelocity = pMomentum / atom_p.mass

    # 设置入射原子的新速度方向
    SetVelocityDirection!(atom_p, pVelocity)

    # 计算入射原子碰撞后的动能
    pEnergy = sum(pMomentum .* pMomentum) / 2 / atom_p.mass

    # 计算总能量转移和局域能量损失
    sumE_t = sum(E_tList)      # 传递给所有目标原子的总动能(单位：eV)
    sumQ_loc = sum(Q_locList)  # 所有局域能量损失总和(单位：eV)

    # 计算可用于动能再分配的能量
    # ENeed = 剩余能量 - 局域能量损失
    ENeed = atom_p.energy - sumQ_loc

    # 计算能量缩放因子λ，确保总能量守恒
    # λ = 可用能量 / (入射原子动能 + 目标原子总动能)
    λ = ENeed / (pEnergy + sumE_t)

    # 将入射原子移动到碰撞点位置
    DisplaceAtom!(atom_p, pPoint, simulator)

    # 设置入射原子缩放后的能量
    SetEnergy!(atom_p, pEnergy * λ)

    # 调试记录（已注释）
    # if atom_p.type == 2
    #     @record "log/$(simulator.nCascade).csv" "$(pEnergy * λ),$(minimum([a.pValue for a in atoms_t])),$(pL),$(N_t),$(atom_p.coordinate[1]),$(atom_p.coordinate[2]),$(atom_p.coordinate[3]),$(atom_p.velocityDirection[1]),$(atom_p.velocityDirection[2]),$(atom_p.velocityDirection[3])" "e,p,pL,N_t,x,y,z,vx,vy,vz" 
    # end

    # 按相同比例缩放所有目标原子的能量
    E_tList *= λ

    # 处理每个目标原子的最终状态
    for (i, atom_t) in enumerate(atoms_t)
        # 检查目标原子是否获得足够能量离开晶格位置
        # 条件：获得能量 > 位移阈值能量 且 净能量 > 0.1 eV
        if E_tList[i] > GetDTE(atom_t, simulator) && E_tList[i] - GetBDE(atom_t, simulator) > 0.1
            # 设置目标原子的净动能（减去结合能）
            SetEnergy!(atom_t, E_tList[i] - GetBDE(atom_t, simulator))
        else
            # 能量不足，目标原子停留在晶格位置，动能为零
            SetEnergy!(atom_t, 0.0)
        end
    end
end

# 动态加载模式下在碰撞级联过程中转储原子状态的函数
function DumpInCascade_dynamicLoad(simulator::Simulator)
    # 检查是否启用了级联过程中的转储功能
    if simulator.parameters.isDumpInCascade
        if simulator.parameters.debugMode == false
            # 生产模式：只转储活动原子和空位
            @dump "Cascade_$(simulator.nCascade).dump" [simulator.atoms; simulator.vacancies] ["vx", "vy", "vz", "e"]
        else
            # 调试模式：转储所有原子，包括晶格原子，用于详细分析
            cells = values(simulator.grid.cells)
            # 收集所有晶格原子
            a = [atom for cell in cells for atom in cell.latticeAtoms]
            # 收集所有活动原子
            b = [atom for cell in cells for atom in cell.atoms]
            # 转储所有原子状态
            @dump "Cascade_$(simulator.nCascade).dump" [a; b] ["vx", "vy", "vz", "e"]
        end
    end
end

# 动态加载模式下的碰撞级联主循环函数
# 模拟入射原子引发的完整碰撞级联过程，包括能量沉积和缺陷产生
function Cascade_dynamicLoad!(atom_p::Atom, simulator::Simulator)
    # 初始化活跃原子列表，开始时只有入射原子
    pAtoms = Vector{Atom}([atom_p])
    # 跟踪所有活跃原子的索引
    pAtomsIndex = [a.index for a in pAtoms]
    parameters = simulator.parameters

    # 初始化碰撞事件计数器
    simulator.nCollisionEvent = 0
    # 增加级联计数器
    simulator.nCascade += 1

    # 转储初始状态
    DumpInCascade_dynamicLoad(simulator)

    # 主碰撞循环：持续直到没有活跃原子
    while true
        # 增加碰撞事件计数
        simulator.nCollisionEvent += 1

        # 初始化存储各原子目标列表的数据结构
        targetsList = Vector{Vector{Atom}}()  # 每个活跃原子的目标原子列表
        deleteIndexes = Int64[]               # 需要删除的原子索引
        othersTargetIndexes = Int64[]         # 其他原子的目标索引（避免重复碰撞）

        # 为每个活跃原子寻找碰撞目标
        for (na, pAtom) in enumerate(pAtoms)
            # 寻找碰撞目标，过滤掉当前活跃原子和已处理的目标
            targets, isAlive = ShotTarget_dynamicLoad(pAtom, [pAtomsIndex; pAtom.lastTargets; othersTargetIndexes], simulator)

            # 检查原子是否存活（未飞出边界）
            if !isAlive
                # 原子飞出边界，清理其目标记录并从系统中删除
                empty!(pAtom.lastTargets)
                delete_dynamicLoad!(simulator, pAtom)
                push!(deleteIndexes, na)  # 标记待删除
                continue
            end

            # 存储找到的目标原子
            push!(targetsList, targets)
            # 记录所有目标原子索引，避免其他原子重复碰撞
            append!(othersTargetIndexes, [t.index for t in targets])
        end

        # 删除飞出边界的原子
        deleteat!(pAtoms, deleteIndexes)
        # 更新活跃原子索引列表
        pAtomsIndex = [a.index for a in pAtoms]

        # 准备下一时间步的活跃原子列表
        nextPAtoms = Vector{Atom}()

        # 处理每个活跃原子的碰撞
        for (pAtom, targets) in zip(pAtoms, targetsList)
            if length(targets) > 0
                # 记录本次碰撞的目标，避免下一时间步重复碰撞
                pAtom.lastTargets = [t.index for t in targets]

                # 执行碰撞计算，更新所有原子的能量和动量
                Collision_dynamicLoad!(pAtom, targets, simulator)

                # 处理获得能量的目标原子
                for target in targets
                    if target.energy > 0.0
                        # 目标原子获得足够能量，离开晶格位置成为反冲原子
                        LeaveLatticePoint_dynamicLoad!(target, simulator)
                        # 位移目标原子到新位置
                        DisplaceAtom!(target, target.coordinate, simulator)
                        # 将反冲原子加入下一时间步的活跃原子列表
                        push!(nextPAtoms, target)
                        # 记录反冲原子的初始碰撞伙伴
                        target.lastTargets = [pAtom.index]
                    end
                end

                # 检查入射原子是否仍有足够能量继续运动
                if pAtom.energy > parameters.stopEnergy
                    # 能量足够，继续参与碰撞级联
                    push!(nextPAtoms, pAtom)
                else
                    # 能量低于停止能量，停止运动并清理目标记录
                    pAtom.lastTargets = Vector{Int64}()
                    Stop_dynamicLoad!(pAtom, simulator)
                end
            else
                # 没有找到碰撞目标，原子继续自由飞行
                push!(nextPAtoms, pAtom)
            end
        end

        # 转储当前时间步的原子状态
        DumpInCascade_dynamicLoad(simulator)

        # 检查是否还有活跃原子继续级联
        if length(nextPAtoms) > 0
            # 更新活跃原子列表，按能量降序排列（优先处理高能原子）
            pAtoms = nextPAtoms
            sort!(pAtoms, by=a -> a.energy, rev=true)
            # 更新活跃原子索引
            pAtomsIndex = [a.index for a in pAtoms]
        else
            # 没有活跃原子，级联结束
            break
        end
    end

    # 内存管理：定期清理晶格原子以减少内存使用
    if simulator.nCascade % parameters.nCascadeEveryLoad == 0
        # 获取当前进程的内存使用量（单位：kB） - 跨平台实现
        # 注释：第二个代码中只有Unix版本，这里保持第一个代码的跨平台实现
        rss = if Sys.iswindows()
            try
                cmd = `wmic process where processid=$(getpid()) get WorkingSetSize /value`
                result = read(cmd, String)
                if occursin("WorkingSetSize=", result)
                    rss_str = split(result, "WorkingSetSize=")[2]
                    rss_str = strip(split(rss_str, "\n")[1])
                    parse(Int, rss_str) ÷ 1024  # 转换为KB
                else
                    0
                end
            catch
                0
            end
        else
            try
                parse(Int, read(`ps -o rss= -p $(getpid())`, String))
            catch
                0
            end
        end
        
        # 检查内存使用是否超过阈值
        if rss > simulator.parameters.maxRSS
            # 执行内存清理，移除未使用的晶格原子
            CleanUpLatticeAtoms(simulator)
        end
    end
end

# 清理晶格原子函数：在动态加载模式下释放内存并重置网格状态
function CleanUpLatticeAtoms(simulator::Simulator)
    # 清空已探索的单元格列表，为下一次碰撞搜索做准备
    empty!(simulator.exploredCells)
    # 清空网格中的所有单元格，释放内存（在动态加载模式下使用字典存储）
    empty!(simulator.grid.cells)
    # 强制进行垃圾回收，释放未使用的内存
    GC.gc()
    # 遍历所有活跃原子，重新分配到对应的网格单元
    for atom in simulator.atoms
        if atom.isAlive
            # 获取原子当前所在的网格单元
            cell = GetCell(simulator.grid, atom.cellIndex)
            # 将原子添加到对应单元的原子列表中
            push!(cell.atoms, atom)
        end
    end
    # 遍历所有活跃空位，重新分配到对应的网格单元
    for vacancy in simulator.vacancies
        if vacancy.isAlive
            # 获取空位当前所在的网格单元
            cell = GetCell(simulator.grid, vacancy.cellIndex)
            # 将空位添加到对应单元的空位列表中
            push!(cell.vacancies, vacancy)
        end
    end
    # 重置最小晶格原子ID计数器为0
    simulator.minLatticeAtomID = 0
end

# 动态加载模式下的原子删除函数：从模拟系统中移除原子
function delete_dynamicLoad!(simulator::Simulator, atom::Atom; isDeleteVacancy::Bool=false)
    # 获取原子所在的网格单元
    cell = GetCell(simulator.grid, atom.cellIndex)
    if !isDeleteVacancy
        # 从单元原子列表中删除该原子
        deleteat!(cell.atoms, findfirst(a -> a.index == atom.index, cell.atoms))
        # 减少模拟系统中的原子总数
        simulator.numberOfAtoms -= 1
    else
        # 从单元空位列表中删除该空位
        deleteat!(cell.vacancies, findfirst(v -> v.index == atom.index, cell.vacancies))
        # 如果单元格已加载晶格原子，需要创建对应的晶格原子
        if cell.isLoaded
            # 将空位类型转换回对应的原子类型（通过减去类型字典长度）
            atom.type -= length(keys(simulator.parameters.typeDict))
            # 创建原子的副本
            latom = CopyAtom(atom, simulator)
            # 分配新的负值ID给晶格原子
            simulator.minLatticeAtomID -= 1
            latom.index = simulator.minLatticeAtomID
            # 将晶格原子添加到单元的晶格原子列表中
            push!(cell.latticeAtoms, latom)
            # 标记为新加载的晶格原子
            latom.isNewlyLoaded = true
        end
        # 减少模拟系统中的空位总数
        simulator.numberOfVacancies -= 1
    end
    # 将原子标记为非活跃状态
    atom.isAlive = false
end

# 动态加载模式下的原子停止函数：处理动能低于停止能量的原子
function Stop_dynamicLoad!(atom::Atom, simulator::Simulator)
    grid = simulator.grid
    # 获取原子所在的网格单元
    cell = GetCell(grid, atom.cellIndex)
    # 初始化最近空位距离平方为无穷大
    nearestVacancyDistance_squared = Inf
    isExist = false
    nearestVacancy = nothing
    # 如果单元的邻近信息未初始化，则进行初始化
    if !cell.isPushedNeighbor
        SetCellNeighborInfo!(cell, grid)
        cell.isPushedNeighbor = true
    end
    # 如果空位恢复距离为0，直接返回（不进行空位恢复）
    if simulator.parameters.vacancyRecoverDistance_squared == 0.0
        return
    end
    # 遍历所有邻近单元格寻找最近的可恢复空位
    for neighborCellInfo in cell.neighborCellsInfo
        index = neighborCellInfo.index
        cross = neighborCellInfo.cross
        # 获取邻近单元格
        neighborCell = GetCell(simulator.grid, index)
        # 遍历邻近单元格中的所有空位
        for vacancy in neighborCell.vacancies
            # 计算原子与空位之间的距离平方（考虑周期性边界条件）
            dr2 = ComputeDistance_squared(atom.coordinate, vacancy.coordinate, cross, simulator.box)
            # 检查空位是否在恢复距离内且是最近的一个
            if dr2 < simulator.parameters.vacancyRecoverDistance_squared && dr2 < nearestVacancyDistance_squared
                nearestVacancyDistance_squared = dr2
                nearestVacancy = vacancy  # 存储最近的空位引用
                isExist = true
            end
        end
    end
    # 如果找到合适的空位，进行空位恢复处理
    if isExist && nearestVacancy !== nothing
        # 检查原子类型是否与空位类型匹配（通过类型偏移量判断）
        if atom.type == nearestVacancy.type - length(keys(simulator.parameters.typeDict))
            # 类型匹配：删除原子和空位（原子填入空位）
            delete_dynamicLoad!(simulator, atom)
            delete_dynamicLoad!(simulator, nearestVacancy, isDeleteVacancy=true)
        else
            # 类型不匹配：将原子移动到空位位置
            SetCoordinate!(atom, nearestVacancy.coordinate)
            # 添加热扰动
            Pertubation!(atom, simulator)
            # 更新原子所在的网格单元
            ChangeCell!(atom, nearestVacancy.cellIndex, simulator)
        end
    end
end

# 动态加载模式下的离开晶格点函数：将原子从晶格位置移出成为间隙原子
function LeaveLatticePoint_dynamicLoad!(atom::Atom, simulator::Simulator; isUpdateEnv::Bool=true)
    # 只对新加载的晶格原子进行处理
    if atom.isNewlyLoaded
        cell = GetCell(simulator.grid, atom.cellIndex)
        # 标记原子已不再是新加载状态
        atom.isNewlyLoaded = false
        # 在原子原来的晶格位置创建空位
        vacancy = CreateVacancy(atom, simulator)
        # 将空位添加到单元的空位列表
        push!(cell.vacancies, vacancy)
        # 将空位添加到模拟器的空位列表
        push!(simulator.vacancies, vacancy)
        # 增加模拟系统中的空位总数
        simulator.numberOfVacancies += 1
        # 为空位分配唯一的ID（使用大数值范围避免与原子ID冲突）
        vacancy.index = simulator.maxVacancyID
        simulator.maxVacancyID += 1

        # 从单元的晶格原子列表中移除该原子
        deleteat!(cell.latticeAtoms, findfirst(a -> a.index == atom.index, cell.latticeAtoms))
        # 将原子添加到单元的原子列表
        push!(cell.atoms, atom)
        # 将原子添加到模拟器的原子列表
        push!(simulator.atoms, atom)
        # 为原子分配新的唯一ID
        simulator.maxAtomID += 1
        atom.index = simulator.maxAtomID
        # 增加模拟系统中的原子总数
        simulator.numberOfAtoms += 1
    end
end

# 原子复制函数：创建原子的浅拷贝副本
function CopyAtom(atom::Atom, simulator::Simulator)
    # 创建新的原子对象，复制类型和坐标
    newAtom = Atom(atom.type, atom.coordinate, simulator.parameters)
    # 复制单元格索引
    newAtom.cellIndex = atom.cellIndex
    return newAtom
end

# 空位创建函数：在原子位置创建对应的空位
function CreateVacancy(atom::Atom, simulator::Simulator)
    # 提取晶格坐标，处理SVector和普通向量的不同情况
    coord = if atom.latticeCoordinate isa SVector
        [atom.latticeCoordinate[1], atom.latticeCoordinate[2], atom.latticeCoordinate[3]]
    else
        atom.latticeCoordinate[:]
    end
    # 创建空位原子对象
    vacancy = Atom(atom.type, coord, simulator.parameters)
    # 复制单元格索引
    vacancy.cellIndex = atom.cellIndex
    # 空位类型 = 原子类型 + 类型字典长度（用于类型区分）
    vacancy.type += length(keys(simulator.parameters.typeDict))
    # 复制在单元格中的索引位置
    vacancy.indexInCell = atom.indexInCell
    return vacancy
end

# 动态加载模式下的目标搜索函数：寻找碰撞目标原子
function ShotTarget_dynamicLoad(atom::Atom, filterIndexes::Vector{Int64}, simulator::Simulator)
    grid = simulator.grid
    periodic = simulator.parameters.periodic
    # 获取原子当前所在的网格单元
    cell = GetCell(grid, atom.cellIndex)
    # 初始化自由飞行路径长度为0
    atom.emptyPath = 0.0
    # 循环搜索直到找到目标或确定无目标
    while true
        # 如果单元的邻近信息未初始化，则进行初始化
        if !cell.isPushedNeighbor
            SetCellNeighborInfo!(cell, grid)
            cell.isPushedNeighbor = true
        end
        # 在当前单元及其邻近单元中搜索碰撞目标
        targets, isInfinity = GetTargetsFromNeighbor_dynamicLoad(atom, cell, filterIndexes, simulator)
        # 如果找到目标原子，返回结果
        if length(targets) > 0
            # 重置所有已探索单元格的标志
            for cell in simulator.exploredCells
                cell.isExplored = false
            end
            # 清空已探索单元格列表
            empty!(simulator.exploredCells)
            return targets, true
        else
            # 没有找到目标，计算原子离开当前单元的面和方向
            dimension, direction, t = AtomOutFaceDimension(atom, cell)
            # 记录自由飞行路径长度
            atom.emptyPath = t
            # 计算邻近单元的索引偏移量（转换为1-based索引）
            neighborIndex = MVector{3,Int8}(0, 0, 0)
            neighborIndex[dimension] = direction == 1 ? Int8(-1) : Int8(1)
            neighborIndex .+= 2
            # 获取邻近单元信息
            neighborInfo = cell.neighborCellsInfo[neighborIndex[1], neighborIndex[2], neighborIndex[3]]
            crossFlag = neighborInfo.cross
            # 处理周期性边界条件：如果穿越边界且该维度是周期性的，调整坐标
            if crossFlag[dimension] != 0 && periodic[dimension]
                atom.coordinate[dimension] -= crossFlag[dimension] * simulator.box.vectors[dimension, dimension]
            end
            # 检查是否到达非周期性边界或无限区域
            if (crossFlag[dimension] != 0 && !periodic[dimension]) || isInfinity
                # 重置已探索单元格状态
                for cell in simulator.exploredCells
                    cell.isExplored = false
                end
                # 重置自由飞行路径长度
                atom.emptyPath = 0.0
                # 清空已探索单元格列表
                empty!(simulator.exploredCells)
                return Vector{Atom}(), false # 返回空结果，表示原子将离开模拟区域
            end
            # 移动到下一个邻近单元继续搜索
            index = neighborInfo.index
            cell = GetCell(grid, index)
        end
    end
end

# 动态加载模式下的数据转储函数：将原子和空位信息写入文件
function Dump_dynamicLoad(simulator::Simulator, fileName::String, step::Int64, type::String="a", isDebug::Bool=false)
    # 检查盒子是否为正交（目前只支持正交盒子的输出）
    if !simulator.parameters.isOrthogonal
        error("The box is not orthogonal, please use the orthogonal box.")
    end
    # 打开文件进行写入
    open(fileName, type) do file
        # 写入LAMMPS数据格式的头部信息
        write(file, "ITEM: TIMESTEP\n")
        write(file, string(step), "\n")
        write(file, "ITEM: NUMBER OF ATOMS\n")
        write(file, string(simulator.numberOfAtoms + simulator.numberOfVacancies), "\n")
        write(file, "ITEM: BOX BOUNDS ")
        # 写入各维度的边界条件类型
        for d in 1:3
            if simulator.parameters.periodic[d]
                write(file, "pp ")  # 周期性边界
            else
                write(file, "ff ")  # 固定边界
            end
        end
        write(file, "\n")
        # 写入各维度的盒子边界
        for d in 1:3
            write(file, "0 $(simulator.box.vectors[d,d])\n")
        end
        # 根据调试模式选择输出格式
        if isDebug
            write(file, "ITEM: ATOMS id type x y z vx vy vz energy cx cy cz dte\n")
        else
            write(file, "ITEM: ATOMS id type x y z e\n")
        end
        # 写入所有活跃原子的信息
        for atom in simulator.atoms
            if atom.isAlive
                if isDebug
                    # 调试模式：输出完整信息包括速度、能量、单元格索引和DTE
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
                    # 普通模式：只输出基本信息和能量
                    write(
                        file,
                        "$(atom.index) $(atom.type) \
            $(atom.coordinate[1]) $(atom.coordinate[2]) $(atom.coordinate[3]) $(atom.energy)\n"
                    )
                end
            end
        end
        # 写入所有活跃空位的信息（ID偏移100000以避免与原子ID冲突）
        for atom in simulator.vacancies
            if atom.isAlive
                if isDebug
                    write(
                        file,
                        "$(atom.index+100000) $(atom.type) \
            $(atom.coordinate[1]) $(atom.coordinate[2]) $(atom.coordinate[3]) \
            $(atom.velocityDirection[1]*sqrt(2*atom.mass*atom.energy)) $(atom.velocityDirection[2]*sqrt(2*atom.mass*atom.energy)) $(atom.velocityDirection[3]*sqrt(2*atom.mass*atom.energy)) \
            $(atom.energy) \
            $(atom.cellIndex[1]) $(atom.cellIndex[2]) $(atom.cellIndex[3]) \
            $(GetDTE(atom, simulator))\n"
                    )
                else
                    write(
                        file,
                        "$(atom.index+100000) $(atom.type) \
            $(atom.coordinate[1]) $(atom.coordinate[2]) $(atom.coordinate[3]) $(atom.energy)\n"
                    )
                end
            end
        end
    end
end

# 动态加载模式下的状态恢复函数：重置模拟器到初始晶格状态
function Restore_dynamicLoad!(simulator::Simulator)
    parameters = simulator.parameters
    # 遍历所有原子和空位
    for atom in [simulator.atoms; simulator.vacancies]
        if atom.isAlive
            cellIndex = atom.cellIndex
            cell = GetCell(simulator.grid, cellIndex)
            # 如果单元格已加载且原子类型是空位类型（通过类型值判断）
            if cell.isLoaded && atom.type > length(keys(simulator.parameters.typeDict))
                # 将空位转换回晶格原子（减去类型偏移量）
                latticeAtom = Atom(atom.type - length(keys(simulator.parameters.typeDict)), atom.coordinate, parameters)
                # 设置晶格坐标为当前坐标
                latticeAtom.latticeCoordinate = SVector{3,Float64}(atom.coordinate[1], atom.coordinate[2], atom.coordinate[3])
                # 添加热扰动
                Pertubation!(latticeAtom, simulator)
                # 设置单元格索引
                latticeAtom.cellIndex = cell.index
                # 分配新的负值ID
                simulator.minLatticeAtomID -= 1
                latticeAtom.index = simulator.minLatticeAtomID
                # 标记为新加载的晶格原子
                latticeAtom.isNewlyLoaded = true
                # 添加到单元的晶格原子列表
                push!(cell.latticeAtoms, latticeAtom)
            end
            # 清空单元的原子列表（在后续步骤中会重新构建）
            empty!(cell.atoms)
        end
    end
    # 清空模拟器的原子和空位列表
    empty!(simulator.atoms)
    empty!(simulator.vacancies)
    # 重置ID计数器
    simulator.maxAtomID = 0
    simulator.maxVacancyID = 1E6
    simulator.numberOfAtoms = 0
    simulator.numberOfVacancies = 0
end