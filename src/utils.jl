# 辅助函数：获取当前线程的随机数生成器
# 优先使用 simulator.workBuffers.threadRNG，如果不可用则回退到 Main.THREAD_RNG（向后兼容）
function get_thread_rng(simulator::Simulator)
    tid = Threads.threadid()
    if hasfield(typeof(simulator.workBuffers), :threadRNG) && 
       length(simulator.workBuffers.threadRNG) >= tid
        return simulator.workBuffers.threadRNG[tid]
    else
        # 向后兼容：如果 WorkBuffers 中没有 threadRNG，使用全局的
        return Main.THREAD_RNG[tid]
    end
end

"""
    DefectStatics(simulator::Simulator)

统计模拟系统中的缺陷（空位和间隙原子）。

返回一个元组 `(interstitials, vacancies)`，其中：
- `interstitials::Vector{Atom}`: 间隙原子列表
- `vacancies::Vector{LatticePoint}`: 空位列表（晶格点对象）

# 参数
- `simulator::Simulator`: 模拟器对象（静态加载模式需要先调用 `Save!`）

# 返回值
返回 `(interstitials, vacancies)` 元组。

# 缺陷识别逻辑
- **空位**: 原始晶格位置没有原子或标记为空位类型
- **间隙原子**: 存活但不在任何晶格点上的原子
- **入射离子**: 存储后添加的原子，如果不在晶格点上也视为间隙原子

# 抛出异常
- `SimulationError`: 如果模拟器未保存（静态加载模式）

# 示例
```julia
# 静态加载模式需要先保存
Save!(simulator)
Cascade!(ion, simulator)
interstitials, vacancies = DefectStatics(simulator)
println("空位数量: ", length(vacancies))
println("间隙原子数量: ", length(interstitials))
```
"""
function DefectStatics(simulator::Simulator)
    # 静态加载模式的缺陷统计
    if !simulator.parameters.is_dynamic_load
        # 检查模拟器是否已保存状态
        if !simulator.isStore
            throw(SimulationError("DefectStatics", "Simulator must be stored when counting defects. Call Save!(simulator) first."))
        end
        
        # 获取晶格点和原子引用
        latticePoints = simulator.latticePoints
        atoms = simulator.atoms
        
        # 初始化缺陷列表
        vacancies = Vector{LatticePoint}()      # 空位列表
        interstitials = Vector{Atom}()          # 间隙原子列表
        
        # 遍历所有位移原子（离开原始位置的原子）
        for idx in simulator.displacedAtoms
            # 检查空位条件：
            # 1. 晶格点没有原子 (atomIndex == -1)
            # 2. 或者晶格点上的原子类型等于类型字典长度（特殊标记，可能是空位类型）
            if latticePoints[idx].atomIndex == -1 || 
               atoms[latticePoints[idx].atomIndex].type == length(keys(simulator.parameters.typeDict))
                push!(vacancies, latticePoints[idx])  # 添加到空位列表
            end     
            # 检查间隙原子条件：
            # 1. 原子存活
            # 2. 不在晶格点上 (latticePointIndex == -1)
            if atoms[idx].isAlive && atoms[idx].latticePointIndex == -1
                push!(interstitials, atoms[idx])  # 添加到间隙原子列表
            end
        end 
        # 遍历存储后添加的原子（主要是入射离子）
        for idx in simulator.numberOfAtomsWhenStored:simulator.maxAtomID
            # 检查间隙原子条件
            if atoms[idx].isAlive && atoms[idx].latticePointIndex == -1
                push!(interstitials, atoms[idx]) # 添加到间隙原子列表
            end
        end  
    else
        # 动态加载模式：直接使用模拟器中的原子和空位列表
        interstitials, vacancies = (simulator.atoms, simulator.vacancies)
    end 
    # 返回缺陷统计结果
    return interstitials, vacancies
end

#简单统计：直接统计所有没有原子的晶格点
function CountVacancies(simulator::Simulator)#空位计数
    nVacancies = 0  # 初始化空位计数器
    # 遍历所有晶格点
    for latticePoint in simulator.latticePoints
        if latticePoint.atomIndex == -1  # 检查晶格点是否为空
            nVacancies += 1  # 空位计数增加
        end
    end 
    return nVacancies  # 返回空位总数
end

#与CountVacancies的区别：返回晶格点对象而非仅计数
function ExtractVacancyLattices(simulator::Simulator)  #提取空位晶格点
    vacancyLattices = []  # 初始化空位晶格点列表 
    # 遍历所有晶格点
    for latticePoint in simulator.latticePoints
        if latticePoint.atomIndex == -1  # 检查是否为空位
            push!(vacancyLattices, latticePoint)  # 添加到空位列表
        end
    end  
    return vacancyLattices  # 返回空位晶格点列表
end

#随机位置生成函数
function RandomPointInCircle(radius::Float64=3.0, simulator::Union{Simulator, Nothing}=nothing)#圆内随机点，均匀分布技巧：sqrt(rand()) 确保点在圆内均匀分布，而非聚集在中心
    rng = if isnothing(simulator)
        Main.THREAD_RNG[Threads.threadid()]  # 向后兼容：如果没有 simulator，使用全局的
    else
        get_thread_rng(simulator)  # 使用 simulator 中的 RNG
    end
    θ = 2π * rand(rng)                    # 生成随机角度 [0, 2π]
    r = sqrt(rand(rng)) * radius          # 生成随机半径，使用sqrt确保均匀分布
    # 转换为笛卡尔坐标
    x = r * cos(θ)
    y = r * sin(θ)
    return [x, y, 0.0]  # 返回二维点（z=0）
end

function RandomInSquare(a::Float64, b::Float64, simulator::Union{Simulator, Nothing}=nothing)#矩形内随机点，应用场景：矩形区域内的均匀采样
    rng = if isnothing(simulator)
        Main.THREAD_RNG[Threads.threadid()]  # 向后兼容：如果没有 simulator，使用全局的
    else
        get_thread_rng(simulator)  # 使用 simulator 中的 RNG
    end
    x = rand(rng) * a  # x坐标在 [0, a] 均匀分布
    y = rand(rng) * b  # y坐标在 [0, b] 均匀分布
    return [x, y, 0.0]  # 返回二维点
end

function RandomInAnUnitGrapheneCell(a::Float64, simulator::Union{Simulator, Nothing}=nothing) #石墨烯原胞内随机点，石墨烯晶格：基于六角晶格的原胞尺寸计算
    # 计算石墨烯原胞尺寸
    X = a * 3          # x方向长度：3倍碳碳键长
    Y = sqrt(3) * a    # y方向长度：√3倍碳碳键长
    rng = if isnothing(simulator)
        Main.THREAD_RNG[Threads.threadid()]  # 向后兼容：如果没有 simulator，使用全局的
    else
        get_thread_rng(simulator)  # 使用 simulator 中的 RNG
    end
    x = rand(rng) * X  # x坐标在 [0, X] 均匀分布
    y = rand(rng) * Y  # y坐标在 [0, Y] 均匀分布
    return [x, y, 0.0] # 二维石墨烯平面
end  


#随机向量生成函数
function RandomVectorInUnitSphere(θrange::Float64, simulator::Union{Simulator, Nothing}=nothing) #单位球内随机向量，θ范围限制：θrange 参数限制向量的最大偏离角度
    rng = if isnothing(simulator)
        Main.THREAD_RNG[Threads.threadid()]  # 向后兼容：如果没有 simulator，使用全局的
    else
        get_thread_rng(simulator)  # 使用 simulator 中的 RNG
    end
    # 球坐标生成
    θ = θrange * rand(rng)  # 极角，限制在 [0, θrange] 范围内
    φ = 2π * rand(rng)      # 方位角，完整 [0, 2π] 范围
    # 转换为笛卡尔坐标
    x = sin(θ) * cos(φ)
    y = sin(θ) * sin(φ)
    z = cos(θ)              # 注意：这里z是正的
    return [x, y, -z]       # 返回向量，z取负（可能指向下方）
end

#=
旋转矩阵计算：
平行向量：单位矩阵
反向向量：180度旋转
一般情况：Rodrigues旋转公式
=#
function rotation_matrix_from_vectors(vec1::AbstractVector, vec2::AbstractVector)#向量间旋转矩阵
    # 归一化输入向量
    a = normalize(vec1)
    b = normalize(vec2)
    # 计算叉积和点积
    v = cross(a, b)        # 旋转轴
    c = dot(a, b)          # 夹角余弦值
    # 特殊情况1：向量平行
    if c ≈ 1.0
        return I(3)        # 返回单位矩阵，无需旋转
    end
    # 特殊情况2：向量反向
    if c ≈ -1.0
        # 选择与a不共线的轴
        other = abs(dot(a, [0, 0, 1])) < 0.9 ? [0, 0, 1] : [1, 0, 0]
        axis = normalize(cross(a, other))  # 计算旋转轴
        
        # 创建180度旋转
        return AngleAxis(π, axis...)|>RotMatrix
    end
    # 一般情况：使用Rodrigues旋转公式
    s = norm(v)                    # 叉积模长
    vx = [0 -v[3] v[2];           # 叉积矩阵
          v[3] 0 -v[1]; 
          -v[2] v[1] 0]
    # Rodrigues公式: R = I + vx + vx²*(1/(1+c))
    R = I(3) + vx + vx^2 * (1 / (1 + c))
    return R
end

#=偏离向量生成流程：
生成高斯分布的极角偏离
生成均匀分布的方位角
在局部坐标系构造偏离向量
计算到目标方向的旋转矩阵
应用旋转得到最终向量=#
function RandomlyDeviatedVector(incident_direction::AbstractVector, divergence::Float64, simulator::Union{Simulator, Nothing}=nothing)#随机偏离向量
    rng = if isnothing(simulator)
        Main.THREAD_RNG[Threads.threadid()]  # 向后兼容：如果没有 simulator，使用全局的
    else
        get_thread_rng(simulator)  # 使用 simulator 中的 RNG
    end
    # 生成随机偏离角度（高斯分布）
    θ = abs(rand(rng, Normal(0.0, divergence)))  # 极角偏离，取绝对值确保非负
    φ = 2π * rand(rng)                           # 随机方位角
    # 在局部坐标系中生成偏离向量
    x = sin(θ) * cos(φ)
    y = sin(θ) * sin(φ)
    z = cos(θ)              # 局部z轴方向
    local_vec = [x, y, -z]  # 局部坐标系中的向量（指向-z）
    # 定义标准轴（指向-z方向）
    standard_axis = [0.0, 0.0, -1.0]
    # 计算从标准轴到入射方向的旋转矩阵
    R = rotation_matrix_from_vectors(standard_axis, incident_direction)
    # 应用旋转，将局部向量转换到入射方向
    final_vec = R * local_vec
    return final_vec
end

# 多层缺陷统计分析函数（慢速实现版本）
# 功能：对模拟系统中的缺陷进行多层结构分析，考虑深度分层和空间分布
# 参数：
#   simulator::Simulator - 模拟器对象，包含所有原子和晶格点信息
#   criteria::Float64 - 缺陷聚类判据，用于确定缺陷聚集的临界距离（单位：Å）
# 返回值：函数体尚未完全实现，框架为未来扩展预留
function MultilayerDefectStatistics_slow(simulator::Simulator, criteria::Float64)
    latticePoints = simulator.latticePoints  # 获取所有晶格点
    vacancies = Vector{LatticePoint}()       # 初始化空位列表
    
    # 第一步：收集所有空位晶格点
    for latticePoint in latticePoints
        if latticePoint.atomIndex == -1  # 判断晶格点是否为空位
            push!(vacancies, latticePoint)  # 添加到空位列表
        end
    end 
    
    # 第二步：分组分析（待实现）
    # 计划功能：
    # 1. 根据深度（z坐标）对缺陷进行分层
    # 2. 每层内进行空间聚类分析
    # 3. 计算各层的缺陷密度和分布
    # 4. 识别缺陷团簇和孤立缺陷
    
    # 函数体待完善...
end

# 基于截断距离和KD树的缺陷聚类分析函数
# 功能：使用KD树空间索引高效地对缺陷坐标进行聚类分析
# 参数：
#   coords::Matrix{Float64} - 缺陷坐标矩阵，尺寸为N×3，每行代表一个缺陷的(x,y,z)坐标
#   rcut::Float64 - 截断距离（单位：Å），用于判断两个缺陷是否属于同一聚类
# 返回值：
#   clusters::Vector{Vector{Int64}} - 聚类列表，每个聚类包含属于该聚类的缺陷索引
function ClusteringByCutoffKDTree(coords::Matrix{Float64}, rcut::Float64)
    N = size(coords, 1)               # 获取缺陷总数
    pts = coords'                     # 转置坐标矩阵，KDTree要求维度×点数格式（3×N）
    tree = KDTree(pts)                # 构建KD树空间索引，加速最近邻搜索
    visited  = falses(N)              # 访问标记数组，初始化为全部未访问
    clusters = Vector{Vector{Int64}}()  # 初始化聚类容器
    
    # 遍历所有缺陷点
    for i in 1:N
        # 如果当前点尚未被访问，则开始一个新的聚类
        if !visited[i]
            queue = [i]                      # 使用队列进行广度优先搜索
            visited[i] = true                # 标记当前点为已访问
            current_cluster = Int64[i]       # 初始化当前聚类，包含起点i
            
            # 广度优先搜索扩展当前聚类
            while !isempty(queue)
                p = pop!(queue)              # 从队列中取出一个点
                # 查找距离点p在rcut范围内的所有邻居点
                neighbors = inrange(tree, pts[:, p], rcut)
                
                # 遍历所有邻居点
                for j in neighbors
                    if !visited[j]           # 如果邻居点尚未访问
                        visited[j] = true    # 标记为已访问
                        push!(queue, j)      # 加入队列，继续扩展搜索
                        push!(current_cluster, j)  # 添加到当前聚类
                    end
                end
            end
            
            # 将完整的聚类添加到聚类列表
            push!(clusters, current_cluster)
        end
    end

    return clusters  # 返回所有聚类
end

# 圆形区域内辐照模拟函数
# 功能：在指定圆形区域内随机生成入射离子，并模拟其碰撞级联过程
# 参数：
#   simulator::Simulator - 模拟器对象
#   energy::Float64 - 入射离子能量（单位：eV）
#   radius::Float64 - 辐照区域半径（单位：Å）
#   center::Vector{Float64} - 辐照区域中心坐标（三维向量，单位：Å）
# 工作流程：
#   1. 在指定圆形区域内随机生成入射离子位置
#   2. 创建入射离子并设置初始速度和能量
#   3. 将离子添加到模拟系统
#   4. 启动碰撞级联模拟
function IrrdiationInCircle(simulator::Simulator, energy::Float64, radius::Float64, center::Vector{Float64})
    # 在指定圆形区域内随机生成入射离子位置
    ionPosition = RandomPointInCircle(radius, simulator) + center  # 使用 simulator 中的 RNG
    
    # 创建入射离子（使用类型字典中的第一种原子类型）
    ion = Atom(simulator.parameters.typeDict[1], ionPosition, simulator.parameters)
    
    # 设置入射离子初始速度方向（垂直向下，沿-z方向）
    SetVelocityDirection!(ion, [0.0, 0.0, -1.0])
    
    # 设置入射离子初始动能
    SetEnergy!(ion, energy)
    
    # 将入射离子添加到模拟系统
    push!(simulator, ion)
    
    # 启动碰撞级联模拟，模拟离子在材料中的能量沉积和缺陷产生
    Cascade!(ion, simulator)
end