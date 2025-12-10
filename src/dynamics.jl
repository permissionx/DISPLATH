# 导入BCA模块和Printf模块，用于二进制碰撞近似计算和格式化输出
using .BCA
using Printf

# 寻找目标原子函数：在模拟器中为入射原子寻找碰撞目标
# 参数：
#   atom::Atom - 入射原子
#   filterIndexes::Vector{Int64} - 需要过滤的原子索引（避免重复碰撞）
#   simulator::Simulator - 模拟器对象
# 返回值：
#   (targets, isAlive, vacancy) - 目标原子向量、原子是否存活标志、可能遇到的空位
function ShotTarget(atom::Atom, filterIndexes::Vector{Int64}, simulator::Simulator)
    grid = simulator.grid  # 获取模拟器的网格系统
    periodic = simulator.parameters.periodic  # 获取各维度的周期性边界条件设置
    cell = GetCell(grid, atom.cellIndex)  # 获取原子当前所在的网格单元
    atom.emptyPath = 0.0  # 初始化自由飞行路径长度为0
    
    # 无限循环，直到找到目标或确定原子飞出系统
    while true
        (targets, isInfinity, vacancy) = GetTargetsFromNeighbor(atom, cell, filterIndexes, simulator)
        
        # 如果找到了空位（且没有找到原子目标）
        if !isnothing(vacancy)
            # 重置所有已探索单元的标记
            for cell in simulator.exploredCells
                cell.isExplored = false
            end
            empty!(simulator.exploredCells)  # 清空已探索单元列表
            return Vector{Atom}(), true, vacancy  # 返回空目标向量、存活标志和找到的空位
        end
        
        # 删除上次碰撞目标中的重复目标（避免重复碰撞） # 新增：优化目标过滤逻辑，避免同一原子在连续时间步中被重复碰撞
        # 如果找到了原子目标
        if length(targets) > 0
            # 重置所有已探索单元的标记
            for cell in simulator.exploredCells
                cell.isExplored = false
            end
            empty!(simulator.exploredCells)  # 清空已探索单元列表
            return targets, true, nothing  # 返回目标原子向量、存活标志和无空位
        else
            # 没有找到目标，计算原子从当前单元射出的面和方向
            dimension, direction, t = AtomOutFaceDimension(atom, cell)
            atom.emptyPath = t  # 记录自由飞行路径长度
            
            # 计算邻近单元的索引偏移量
            neighborIndex = Vector{Int8}([0,0,0])
            neighborIndex[dimension] = direction == 1 ? Int8(-1) : Int8(1)  # 根据射出方向确定偏移
            neighborIndex .+= 2  # 将索引转换为1-based（从-1,0,1转换为1,2,3）
            
            # 获取邻近单元信息
            neighborInfo = cell.neighborCellsInfo[neighborIndex...]
            crossFlag = neighborInfo.cross  # 获取边界穿越标志
            
            # 如果存在周期性边界穿越且该维度是周期性的
            if crossFlag[dimension] != 0 && periodic[dimension]
                # 应用周期性边界条件：将原子坐标调整到对边
                atom.coordinate[dimension] -= crossFlag[dimension] * simulator.box.vectors[dimension, dimension]
            end
            
            # 如果遇到非周期性边界或无限飞行情况
            if (crossFlag[dimension] != 0 && !periodic[dimension]) || isInfinity
                # 重置所有已探索单元的标记
                for cell in simulator.exploredCells
                    cell.isExplored = false
                end
                atom.emptyPath = 0.0  # 重置自由飞行路径长度
                empty!(simulator.exploredCells)  # 清空已探索单元列表
                return Vector{Atom}(), false, nothing  # 返回无目标、原子不存活、无空位
            end 
            
            # 移动到下一个邻近单元继续搜索
            index = neighborInfo.index
            cell = GetCell(grid, index)
        end
    end
end

# 计算原子从网格单元射出的维度和方向
# 参数：
#   atom::Atom - 入射原子
#   cell::Cell - 当前网格单元
# 返回值：
#   (dimension, direction, t) - 射出维度、方向(1表示min面，2表示max面)、飞行时间
function AtomOutFaceDimension(atom::Atom, cell::Cell)
    coordinate = atom.coordinate  # 获取原子当前位置
    
    # 遍历三个维度
    for d in 1:3
        # 根据速度方向确定要检查的单元面
        if atom.velocityDirection[d] >= 0
            rangeIndex = 2  # 正向速度检查max面
        else
            rangeIndex = 1  # 负向速度检查min面
        end
        
        faceCoordinate = cell.ranges[d, rangeIndex]  # 获取单元面的坐标位置
        
        # 计算到达该面的飞行时间：t = (面坐标 - 原子坐标) / 速度分量
        t = (faceCoordinate - coordinate[d]) / atom.velocityDirection[d]
        
        # 获取其他两个维度的索引
        elseDs = [ed for ed in 1:3 if ed != d]
        allInRange = true  # 标记其他维度坐标是否在单元范围内
        
        # 检查在其他维度上，飞行时间t后坐标是否仍在单元内
        for elseD in elseDs
            crossCoord = coordinate[elseD] + atom.velocityDirection[elseD] * t  # 计算交叉点坐标
            if !(cell.ranges[elseD, 1] <= crossCoord <= cell.ranges[elseD, 2])
                allInRange = false  # 如果超出范围，标记为false
                break
            end
        end
        
        # 如果所有维度的交叉点都在单元内，返回当前维度信息
        if allInRange
            return d, rangeIndex, t
        end
    end
    
    # 如果找不到射出面，抛出错误（理论上不应该发生）
    error("Out face not found\n ########Atom#######\n $(atom) \n 
                                ########cell#######\n $(cell.ranges) \n $(cell.index)\n")
end

# 从邻近单元中获取可能的碰撞目标
# 参数：
#   atom::Atom - 入射原子
#   cell::Cell - 当前网格单元
#   filterIndexes::Vector{Int64} - 需要过滤的原子索引
#   simulator::Simulator - 模拟器对象
# 返回值：
#   (targets, infiniteFlag, nearestVacancy) - 目标原子、无限飞行标志、最近空位
function GetTargetsFromNeighbor(atom::Atom, cell::Cell, filterIndexes::Vector{Int64}, simulator::Simulator)
    grid = simulator.grid  # 获取网格系统
    box = simulator.box    # 获取模拟盒子
    targets = Vector{Atom}()  # 初始化目标原子向量
    infiniteFlag = true    # 初始化无限飞行标志（假设没有找到任何单元）
    candidateTargets = Vector{Atom}()  # 初始化候选目标原子向量
    pMax = simulator.parameters.pMax  # 获取最大碰撞参数
    minVacancyPL = Float64(Inf)  # 初始化最小空位路径长度为无穷大
    nearestVacancy::Union{Atom, Nothing} = nothing  # 初始化最近空位为nothing
    
    # 遍历所有邻近单元（3x3x3=27个方向）
    for neighborCellInfo in cell.neighborCellsInfo
        cross = neighborCellInfo.cross  # 获取边界穿越信息
        nonPeriodicFlag = false  # 初始化非周期性边界标志
        
        # 检查是否存在非周期性边界的穿越
        for d in 1:3
            if cross[d] != 0 && !simulator.parameters.periodic[d] 
                nonPeriodicFlag = true  # 标记存在非周期性边界穿越
                break
            end
        end
        
        # 如果存在非周期性边界穿越，跳过该邻近单元
        if nonPeriodicFlag
            continue 
        end 
        
        index = neighborCellInfo.index  # 获取邻近单元索引
        neighborCell = GetCell(grid, index)  # 获取邻近单元对象
        
        # 如果该单元已被探索过，跳过
        if neighborCell.isExplored
            continue
        end
        
        neighborCell.isExplored = true  # 标记该单元为已探索
        push!(simulator.exploredCells, neighborCell)  # 添加到已探索单元列表
        infiniteFlag = false  # 找到至少一个单元，清除无限飞行标志
        
        # 遍历邻近单元中的所有原子
        for neighborAtom in neighborCell.atoms
            # 跳过自身和需要过滤的原子
            if neighborAtom.index == atom.index || neighborAtom.index in filterIndexes    
                continue
            end
            
            # 计算原子与目标原子在速度方向上的距离（正数表示在前方）
            if ComputeVDistance(atom, neighborAtom, neighborCellInfo.cross, box) > 0 
                p = ComputeP!(atom, neighborAtom, neighborCellInfo.cross, box)  # 计算碰撞参数
                
                # 如果碰撞参数超过最大值，跳过该原子
                if p >= pMax
                    continue
                end
                
                push!(candidateTargets, neighborAtom)  # 添加到候选目标
            end
        end
        
        # 如果原子能量低于位移阈值能量，检查空位碰撞
        if atom.energy <= GetDTE(atom, simulator)
            # 遍历邻近单元中的所有空位
            for vacancy in neighborCell.vacancies
                # 计算原子与空位在速度方向上的距离
                if ComputeVDistance(atom, vacancy, neighborCellInfo.cross, box) > 0 
                    p = ComputeP!(atom, vacancy, neighborCellInfo.cross, box)  # 计算碰撞参数
                    
                    # 如果碰撞参数超过最大值，跳过该空位
                    if p >= pMax  # 需要在参数中分配此值
                        continue
                    end 
                    
                    # 更新最近空位信息
                    if vacancy.pL <= minVacancyPL
                        minVacancyPL = vacancy.pL
                        nearestVacancy = vacancy
                    end
                end
            end
        end
    end
    
    # 如果标记为无限飞行（未找到任何有效单元），记录日志
    if infiniteFlag
        @record "log" "Infinitely fly atom in the $(simulator.nCascade)th irradiation:\n$(atom)"
    end  # 无限检查应该在ShotTarget函数中应用
    
    # 如果没有找到候选目标原子
    if isempty(candidateTargets)
        return targets, infiniteFlag, nearestVacancy
    end
    
    # 在候选目标中找到路径长度最小的目标（最近的碰撞目标）
    _, minIdx = findmin(t -> t.pL, candidateTargets)
    nearestTarget = candidateTargets[minIdx]  
    
    # 如果存在空位且空位比原子目标更近，返回空位信息
    if !isnothing(nearestVacancy) && nearestVacancy.pL <= nearestTarget.pL
        return targets, infiniteFlag, nearestVacancy 
    end
    
    push!(targets, nearestTarget)  # 添加最近目标到结果向量
    
    # 检查其他候选目标是否满足同时碰撞条件
    for candidateTarget in candidateTargets
        if candidateTarget.index == nearestTarget.index
            continue  # 跳过自身
        end
        
        matchFlag = true  # 初始化匹配标志
        # 检查与所有已选目标的兼容性
        for target in targets            
            if !SimultaneousCriteria(candidateTarget, target, simulator)
                matchFlag = false  # 不满足同时碰撞条件
                break
            end
        end
        
        # 如果满足条件，添加到目标列表
        if matchFlag
            push!(targets, candidateTarget)
        end
    end    
    
    return (targets, infiniteFlag, nothing)  # 返回目标原子、无限飞行标志、无空位
end

# 碰撞处理函数：计算入射原子与多个目标原子的碰撞
# 参数：
#   atom_p::Atom - 入射原子
#   atoms_t::Vector{Atom} - 目标原子向量
#   simulator::Simulator - 模拟器对象
function Collision_!(atom_p::Atom, atoms_t::Vector{Atom}, simulator::Simulator)
    N_t = length(atoms_t)  # 目标原子数量
    grid = simulator.grid  # 获取网格系统
    
    # 初始化碰撞参数列表
    tanφList = Vector{Float64}(undef, N_t)  # 入射原子散射角正切值
    tanψList = Vector{Float64}(undef, N_t)  # 目标原子散射角正切值
    E_tList = Vector{Float64}(undef, N_t)   # 传递给目标原子的能量
    x_pList = Vector{Float64}(undef, N_t)   # 入射原子位移
    x_tList = Vector{Float64}(undef, N_t)   # 目标原子位移
    Q_locList = Vector{Float64}(undef, N_t) # 局域能量损失
    
    atom_t = atoms_t[1]  # 获取第一个目标原子（用于参考计算）
    pL = atom_t.pL - atom_p.emptyPath  # 计算有效碰撞路径长度
    pPoint = atom_t.pPoint  # 获取碰撞点坐标
    N = GetCell(grid, atom_t.cellIndex).atomicDensity  # 获取原子数密度
    
    # 计算非局域能量损失（如果启用）
    if !simulator.parameters.isNonQnl
        Q_nl_v = Q_nl(atom_p.energy, atom_p.mass, atom_t.mass, atom_p.type, atom_t.type,
                            pL, N, simulator.constantsByType)  
        atom_p.energy -= Q_nl_v  # 从入射原子能量中减去非局域损失
    else
        Q_nl_v = 0.0  # 禁用非局域能量损失
    end
    
    # 能量修正：确保原子能量不会略低于停止能量
    if atom_p.energy < 0.1 && atom_p.energy + Q_nl_v >= 0.1
        atom_p.energy = 0.11
    end
    
    momentum = Vector{Float64}([0.0,0.0,0.0])  # 初始化动量向量
    
    # 遍历所有目标原子，计算碰撞参数
    for (i, atom_t) in enumerate(atoms_t) 
        p = atom_t.pValue  # 获取碰撞参数值
        
        # 计算碰撞物理参数
        tanφList[i], tanψList[i], E_tList[i], x_pList[i], x_tList[i], Q_locList[i] = CollisionParams(
            atom_p.energy, atom_p.mass, atom_t.mass, atom_p.type, atom_t.type, p, simulator.constantsByType,
            simulator.θFunctions[[atom_p.type, atom_t.type]], simulator.τFunctions[[atom_p.type, atom_t.type]])
        
        # 计算目标原子的新速度方向
        if atom_t.pValue != 0
            velocityDirectionTmp = -atom_t.pVector / atom_t.pValue * tanψList[i] + atom_p.velocityDirection
        else
            velocityDirectionTmp = atom_p.velocityDirection  # 零碰撞参数情况
        end   
        
        SetVelocityDirection!(atom_t, velocityDirectionTmp)  # 设置目标原子速度方向
        momentum += sqrt(2 * atom_t.mass * E_tList[i]) * atom_t.velocityDirection  # 累加动量
    end
    
    # 计算入射原子的新动量和速度
    pMomentum = sqrt(2 * atom_p.mass * atom_p.energy) * atom_p.velocityDirection - momentum
    pVelocity = pMomentum  / atom_p.mass
    SetVelocityDirection!(atom_p, pVelocity)  # 设置入射原子新速度方向
    
    # 计算入射原子的剩余能量
    pEnergy =  sum(pMomentum .* pMomentum) / 2 / atom_p.mass
    sumE_t = sum(E_tList)      # 总传递给目标原子的能量
    sumQ_loc = sum(Q_locList)  # 总局域能量损失
    
    # 能量守恒计算
    ENeed = atom_p.energy - sumQ_loc  # 入射原子可用于动能分配的能量
    λ = ENeed / (pEnergy + sumE_t)    # 能量分配比例因子
    
    # 位移入射原子到碰撞点
    DisplaceAtom!(atom_p, pPoint, simulator)
    SetEnergy!(atom_p, pEnergy * λ)  # 设置入射原子新能量
    
    # 按比例调整目标原子获得的能量
    E_tList *= λ
    
    # 处理每个目标原子的状态
    for (i, atom_t) in enumerate(atoms_t)
        # 如果目标原子获得的能量超过其位移阈值和结合能
        if E_tList[i] > GetDTE(atom_t, simulator) && E_tList[i] - GetBDE(atom_t, simulator) > 0.1
            SetEnergy!(atom_t, E_tList[i] - GetBDE(atom_t, simulator))  # 设置目标原子能量
            atom_t.pAtomIndex = atom_p.index  # 临时存储：记录入射原子索引
            atom_t.pDirection = atom_p.velocityDirection  # 临时存储：记录入射方向
        else
            # 能量不足，目标原子停止运动
            SetEnergy!(atom_t, 0.0)
            SetVelocityDirection!(atom_t, SVector{3,Float64}([0.0,0.0,0.0]))
        end
    end
end 

# 碰撞处理函数：处理入射原子与多个目标原子之间的碰撞动力学
# 输入参数：
#   atom_p::Atom - 入射原子（初级碰撞原子）
#   atoms_t::Vector{Atom} - 目标原子向量（可能同时与多个原子碰撞）
#   simulator::Simulator - 模拟器对象，包含所有模拟状态和参数
function Collision!(atom_p::Atom, atoms_t::Vector{Atom}, simulator::Simulator)
    N_t = length(atoms_t)  # 目标原子数量
    grid = simulator.grid  # 模拟网格对象
    
    # 预分配碰撞参数数组，存储每个目标原子的碰撞计算结果
    tanφList = Vector{Float64}(undef, N_t)  # 入射原子散射角正切值数组
    tanψList = Vector{Float64}(undef, N_t)  # 目标原子散射角正切值数组
    E_tList = Vector{Float64}(undef, N_t)   # 传递给每个目标原子的能量数组(单位：eV)
    x_pList = Vector{Float64}(undef, N_t)   # 入射原子位移距离数组(单位：Å)
    x_tList = Vector{Float64}(undef, N_t)   # 目标原子位移距离数组(单位：Å)
    Q_locList = Vector{Float64}(undef, N_t) # 局域能量损失数组(单位：eV)
    
    atom_t = atoms_t[1]  # 取第一个目标原子作为参考（用于密度计算）
    pL = atom_t.pL - atom_p.emptyPath  # 有效碰撞路径长度 = 几何路径 - 自由飞行路径(单位：Å)
    pPoint = atom_t.pPoint  # 碰撞点空间坐标(单位：Å)
    N = GetCell(grid, atom_t.cellIndex).atomicDensity  # 目标原子所在网格单元的原子数密度(单位：原子/Å³)
    
    # 计算非局域能量损失（电子阻止本领），除非明确关闭
    if !simulator.parameters.isNonQnl
        # 调用Q_nl函数计算非局域能量损失
        Q_nl_v = Q_nl(atom_p.energy, atom_p.mass, atom_t.mass, atom_p.type, atom_t.type,
                            pL, N, simulator.constantsByType)  
        atom_p.energy -= Q_nl_v  # 从入射原子能量中扣除非局域能量损失
    else
        Q_nl_v = 0.0  # 如果关闭非局域能量损失，设为0
    end
    
    # 能量修正：如果扣除能量损失后能量低于阈值但原本高于阈值，设为略高于阈值
    if atom_p.energy < 0.1 && atom_p.energy + Q_nl_v >= 0.1
        atom_p.energy = 0.11  # 设为0.11 eV，略高于停止能量阈值
    end
    
    # 遍历所有目标原子，计算每个原子的碰撞参数
    for (i, atom_t) in enumerate(atoms_t) 
        p = atom_t.pValue  # 当前目标原子的碰撞参数p(单位：Å)
        #N = GetCell(grid, atom_t.cellIndex).atomicDensity   # 注释：可选的局部密度计算
        
        # 调用CollisionParams函数计算完整的碰撞参数
        tanφList[i], tanψList[i], E_tList[i], x_pList[i], x_tList[i], Q_locList[i] = CollisionParams(
            atom_p.energy, atom_p.mass, atom_t.mass, atom_p.type, atom_t.type, p, simulator.constantsByType,
            simulator.θFunctions[[atom_p.type, atom_t.type]], simulator.τFunctions[[atom_p.type, atom_t.type]])
    end
    
    sumE_t = sum(E_tList)     # 传递给所有目标原子的总能量(单位：eV)
    sumQ_loc = sum(Q_locList) # 所有目标原子的总局域能量损失(单位：eV)
    
    # 计算能量分配因子η，确保能量守恒
    # η = (入射原子能量) / (入射原子能量 + 其他目标原子的能量损失)
    η = N_t * atom_p.energy / (N_t * atom_p.energy + (N_t - 1) * (sumE_t + sumQ_loc))
    E_tList *= η  # 按比例缩放传递给目标原子的能量
    
    #avePPoint = Vector{Float64}([0.0,0.0,0.0])  # 注释：平均碰撞点坐标计算
    momentum = Vector{Float64}([0.0,0.0,0.0])  # 动量向量，用于动量守恒计算

    # 处理每个目标原子的碰撞后状态
    for (i, atom_t) in enumerate(atoms_t)
        # 计算目标原子的新速度方向
        if atom_t.pValue != 0
            # 根据碰撞几何计算新的速度方向
            velocityDirectionTmp = -atom_t.pVector / atom_t.pValue * tanψList[i] + atom_p.velocityDirection
        else
            # 如果碰撞参数为0（对心碰撞），保持入射原子方向
            velocityDirectionTmp = atom_p.velocityDirection
        end   
        SetVelocityDirection!(atom_t, velocityDirectionTmp)  # 设置目标原子速度方向
        
        #avePPoint += atom_t.pPoint  # 注释：累加碰撞点坐标
        
        # 累加目标原子的动量贡献
        momentum += sqrt(2 * atom_t.mass * E_tList[i]) * atom_t.velocityDirection
        
        # 判断目标原子是否被撞离晶格位置
        if E_tList[i] > GetDTE(atom_t, simulator) && E_tList[i] - GetBDE(atom_t, simulator) > 0.1
            # 如果传递能量超过位移阈值和结合能，目标原子获得动能
            SetEnergy!(atom_t, E_tList[i] - GetBDE(atom_t, simulator))
            
            #tCoordinate = atom_t.coordinate + x_tList[i] * η * atom_p.velocityDirection  # 注释：目标原子位移计算
            #DisplaceAtom!(atom_t, tCoordinate, simulator)  # 注释：目标原子位置更新
            
            # 临时存储关联信息，用于后续处理
            atom_t.pAtomIndex = atom_p.index        # 存储关联的入射原子索引
            atom_t.pDirection = atom_p.velocityDirection  # 存储入射方向
        else 
            # 如果能量不足，目标原子停止运动
            SetEnergy!(atom_t, 0.0)  # 能量设为0
            SetVelocityDirection!(atom_t, SVector{3,Float64}([0.0,0.0,0.0]))  # 速度方向归零
        end
    end

    # 计算入射原子的新位置和状态
    #avePPoint /= N_t  # 注释：计算平均碰撞点
    #x_p = η * sum(x_pList) / N_t  # 注释：计算入射原子位移距离
    #pCoordinate = avePPoint - x_p * atom_p.velocityDirection  # 注释：计算入射原子新位置
    
    # 将入射原子移动到碰撞点位置
    DisplaceAtom!(atom_p, pPoint, simulator)
    
    # 根据动量守恒计算入射原子的新速度
    # 新速度 = (入射原子原动量 - 目标原子总动量) / 入射原子质量
    velocity = (sqrt(2 * atom_p.mass * atom_p.energy) * atom_p.velocityDirection - momentum) / atom_p.mass
    SetVelocityDirection!(atom_p, velocity)  # 设置入射原子新速度方向
    
    # 更新入射原子能量：扣除传递给目标原子的能量和局域能量损失
    SetEnergy!(atom_p, atom_p.energy - (sumE_t + sumQ_loc) * η)
end 

# 碰撞级联主函数：根据加载模式选择相应的级联模拟方法
# 输入参数：
#   atom_p::Atom - 入射原子
#   simulator::Simulator - 模拟器对象
function Cascade!(atom_p::Atom, simulator::Simulator)
    if !IS_DYNAMIC_LOAD
        Cascade_staticLoad!(atom_p, simulator)  # 静态加载模式
    else
        Cascade_dynamicLoad!(atom_p, simulator)  # 动态加载模式
    end
end

# 静态加载模式下的碰撞级联模拟
# 输入参数：
#   atom_p::Atom - 入射原子
#   simulator::Simulator - 模拟器对象
function Cascade_staticLoad!(atom_p::Atom, simulator::Simulator)
    pAtoms = Vector{Atom}([atom_p])  # 初始化活动原子列表，开始时只有入射原子
    pAtomsIndex = [a.index for a in pAtoms]  # 活动原子的索引列表
    parameters = simulator.parameters  # 模拟参数
    
    simulator.nCollisionEvent = 0  # 重置碰撞事件计数器
    simulator.nCascade += 1        # 增加级联计数器
    
    DumpInCascade(simulator)  # 输出初始状态（如果启用）
    
    # 主循环：持续进行直到没有活动原子
    while true
        simulator.nCollisionEvent += 1  # 增加碰撞事件计数
        
        targetsList = Vector{Vector{Atom}}()     # 存储每个活动原子的目标原子列表
        deleteIndexes = Int64[]                  # 需要删除的活动原子索引
        othersTargetIndexes = Int64[]            # 其他活动原子的目标索引（避免重复碰撞）
        
        # 遍历所有活动原子，寻找碰撞目标
        for (na, pAtom) in enumerate(pAtoms)
            # 寻找碰撞目标，过滤掉自身和最近碰撞过的原子
            targets, isAlive, vacancy = ShotTarget(pAtom, [pAtomsIndex; pAtom.lastTargets; othersTargetIndexes], simulator)
            
            # 检查是否遇到空位（vacancy recovery）
            if !isnothing(vacancy) 
                # 如果遇到空位，将活动原子放置在空位位置
                latticePoint = simulator.latticePoints[vacancy.index]
                SetOnLatticePoint!(pAtom, latticePoint, simulator)
                
                # 从系统中删除该空位
                deleteat!(simulator.vacancies, findfirst(v -> v.index == vacancy.index, simulator.vacancies))
                cell = GetCell(simulator.grid, latticePoint.cellIndex)
                deleteat!(cell.vacancies, findfirst(v -> v.index == vacancy.index, cell.vacancies))
                
                empty!(pAtom.lastTargets)  # 清空上次碰撞目标记录
                push!(deleteIndexes, na)   # 标记该原子需要从活动列表删除
                continue 
            end
            
            # 检查原子是否存活（是否飞出边界等）
            if !isAlive
                empty!(pAtom.lastTargets)  # 清空碰撞记录
                delete!(simulator, pAtom)  # 从模拟器中删除原子
                push!(deleteIndexes, na)   # 标记删除
                continue
            end
            
            push!(targetsList, targets)  # 保存该原子的目标列表
            append!(othersTargetIndexes, [t.index for t in targets])  # 记录所有目标索引
        end
        
        # 删除已停止或遇到空位的原子
        deleteat!(pAtoms, deleteIndexes)
        pAtomsIndex = [a.index for a in pAtoms]  # 更新活动原子索引列表
        
        nextPAtoms = Vector{Atom}()  # 下一时间步的活动原子列表
        
        # 处理每个活动原子与其目标的碰撞
        for (pAtom, targets) in zip(pAtoms, targetsList)
            if length(targets) > 0  # 如果有碰撞目标
                pAtom.lastTargets = [t.index for t in targets]  # 记录本次碰撞的目标
                Collision!(pAtom, targets, simulator)  # 执行碰撞计算
                
                # 处理碰撞后的目标原子
                for target in targets
                    if target.energy > 0.0   # 如果目标原子获得足够能量
                        # 更新目标原子位置以确保其在正确的网格单元中，即使坐标未改变也需要调用此函数来更新cell索引
                        DisplaceAtom!(target, target.coordinate, simulator)  # 新增：确保目标原子在正确的网格单元中
                        push!(nextPAtoms, target)  # 将目标原子加入下一时间步活动列表
                        target.lastTargets = [pAtom.index]  # 记录碰撞来源
                        
                        # 如果目标原子原本在晶格位置上，将其移出晶格
                        if target.latticePointIndex != -1
                            LeaveLatticePoint!(target, simulator)
                        end    
                    end
                end
                
                # 处理碰撞后的入射原子
                if pAtom.energy > parameters.stopEnergy 
                    push!(nextPAtoms, pAtom)  # 如果仍有能量，继续运动
                else
                    pAtom.lastTargets = Vector{Int64}()  # 清空碰撞记录
                    Stop!(pAtom, simulator)  # 停止原子运动
                end
            else      
                push!(nextPAtoms, pAtom)  # 如果没有碰撞目标，原子继续直线运动
            end
        end
        
        DumpInCascade(simulator)  # 输出当前状态（如果启用）
        
        # 检查是否还有活动原子
        if length(nextPAtoms) > 0
            pAtoms = nextPAtoms  # 更新活动原子列表
            # 按能量降序排序，优先处理高能原子
            sort!(pAtoms, by = a -> a.energy, rev = true)
            pAtomsIndex = [a.index for a in pAtoms]  # 更新索引列表
        else
            break  # 没有活动原子，结束级联
        end
    end
end

# 级联过程中输出原子状态函数
# 输入参数：
#   simulator::Simulator - 模拟器对象
function DumpInCascade(simulator::Simulator)
    if simulator.parameters.isDumpInCascade  # 检查是否启用级联输出
        # 使用@dump宏输出原子状态到文件
        # 文件名格式：Cascade_{级联编号}.dump
        # 输出属性：速度分量(vx,vy,vz)和能量(e)
        @dump "Cascade_$(simulator.nCascade).dump" simulator.atoms ["vx", "vy", "vz", "e"]
    end
end