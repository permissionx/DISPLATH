# 导入必要的Julia包
#using PyCall  # 注释掉的Python调用接口，可能用于未来扩展

"""
    Box

表示模拟系统的空间边界和周期性条件。

# 字段
- `vectors::Matrix{Float64}`: 3x3矩阵，盒子的基向量(单位：Å)，定义模拟空间的三个方向
- `reciprocalVectors::Matrix{Float64}`: 3x3矩阵，倒易空间基向量(单位：Å⁻¹)，用于周期性边界条件计算
- `isOrthogonal::Bool`: 布尔标志，指示盒子是否为正交坐标系（基向量相互垂直）

# 构造函数
```julia
Box(vectors::Matrix{Float64})
```

# 示例
```julia
# 创建一个正交盒子（10×10×10 Å）
vectors = [10.0 0.0 0.0; 0.0 10.0 0.0; 0.0 0.0 10.0]
box = Box(vectors)
```
"""
mutable struct Box
    vectors::Matrix{Float64}          # 3x3矩阵，表示盒子的基向量(单位：Å)
    reciprocalVectors::Matrix{Float64} # 3x3矩阵，表示倒易空间基向量(单位：Å⁻¹)
    isOrthogonal::Bool                # 布尔标志，指示盒子是否为正交坐标系
end

"""
    Atom

表示模拟系统中的单个原子。

# 字段
- `index::Int64`: 原子的唯一标识符，在整个模拟过程中保持不变
- `isAlive::Bool`: 布尔标志，指示原子是否处于活动状态（是否参与模拟）
- `type::Int64`: 原子类型标识符，对应typeDict中的键值
- `coordinate::Vector{Float64}`: 3维向量，原子的当前位置坐标(单位：Å)
- `cellIndex::Tuple{Int64, Int64, Int64}`: 三元组，表示原子所在网格单元的索引(x,y,z)
- `radius::Float64`: 原子的有效半径(单位：Å)，用于碰撞检测
- `mass::Float64`: 原子质量(单位：原子质量单位u)
- `velocityDirection::SVector{3,Float64}`: 3维静态向量，速度方向单位向量(无量纲)
- `energy::Float64`: 原子的当前动能(单位：eV)
- `Z::Float64`: 原子序数，表示核电荷数
- `dte::Float64`: 位移阈值能量(Displacement Threshold Energy)，单位：eV
- `bde::Float64`: 结合能(Binding Energy)，单位：eV
- `emptyPath::Float64`: 自由飞行路径长度，在碰撞前穿越的距离(单位：Å)

# 构造函数
```julia
Atom(type::Int64, coordinate::Vector{Float64}, parameters::Parameters)
```

# 示例
```julia
# 创建一个碳原子在原点
atom = Atom(1, [0.0, 0.0, 0.0], parameters)
SetEnergy!(atom, 100.0)  # 设置动能为 100 eV
SetVelocityDirection!(atom, [0.0, 0.0, -1.0])  # 设置速度方向向下
```
"""
mutable struct Atom
    index::Int64                      # 原子的唯一标识符，在整个模拟过程中保持不变
    isAlive::Bool                     # 布尔标志，指示原子是否处于活动状态（是否参与模拟）
    type::Int64                       # 原子类型标识符，对应typeDict中的键值
    coordinate::Vector{Float64}       # 3维向量，原子的当前位置坐标(单位：Å)
    cellIndex::Tuple{Int64, Int64, Int64}  # 三元组，表示原子所在网格单元的索引(x,y,z)
    radius::Float64                   # 原子的有效半径(单位：Å)，用于碰撞检测
    mass::Float64                     # 原子质量(单位：原子质量单位u)
    velocityDirection::SVector{3,Float64}  # 3维静态向量，速度方向单位向量(无量纲)
    energy::Float64                   # 原子的当前动能(单位：eV)
    Z::Float64                        # 原子序数，表示核电荷数

    dte::Float64                      # 位移阈值能量(Displacement Threshold Energy)，单位：eV
    bde::Float64                      # 结合能(Binding Energy)，单位：eV

    #numberOfEmptyCells::Int64        # 注释：穿越的空单元数量，用于路径长度计算
    emptyPath::Float64                # 自由飞行路径长度，在碰撞前穿越的距离(单位：Å)

    # 以下字段用于目标原子(atom_t)的碰撞计算
    pValue::Float64                   # 碰撞参数p的值，即碰撞点到目标原子连线在垂直方向的距离(单位：Å)
    pPoint::SVector{3,Float64}        # 3维静态向量，碰撞点的空间坐标(单位：Å)
    pVector::SVector{3,Float64}       # 3维静态向量，从目标原子指向碰撞点的向量(单位：Å)
    pL::Float64                       # 碰撞路径长度参数，用于非局域能量损失计算(单位：Å)
    pAtomIndex::Int64                 # 临时字段，关联的入射原子索引
    pDirection::Vector{Float64}       # 临时字段，碰撞方向向量

    # 以下字段用于入射原子(atom_p)的碰撞计算
    lastTargets::Vector{Int64}        # 整数向量，存储上一次碰撞中涉及的目标原子索引

    latticePointIndex::Int64          # 晶格点索引，-1表示原子不在晶格位置上
    
    # 以下字段用于动力学蒙特卡洛(KMC)模拟
    frequency::Float64                # 迁移频率(单位：Hz或s⁻¹)
    frequencies::Vector{Float64}      # 浮点数向量，存储向不同晶格点迁移的频率
    finalLatticePointIndexs::Vector{Int64}  # 整数向量，可迁移到的目标晶格点索引
    eventIndex::Int64                 # 在KMC事件列表中的索引

    # 以下字段用于动态加载系统
    isNewlyLoaded::Bool               # 布尔标志，指示原子是否为新加载的晶格原子
    latticeCoordinate::SVector{3,Float64}  # 3维静态向量，原子的理想晶格位置坐标(单位：Å)
    indexInCell::Int64                # 原子在所在网格单元中的局部索引
end

"""
    Material

表示完整的材料系统，包含几何结构、原子信息和空间网格定义。

该结构体提供了对模拟系统的高层抽象，便于封装和传递材料数据。

# 字段
- `box::Box`: 模拟盒子对象，定义材料的空间边界和周期性条件
- `atoms::Vector{Atom}`: 原子向量，存储材料中的所有原子（静态加载模式）
- `inputGridVectors::Matrix{Float64}`: 3x3矩阵，输入网格向量，用于定义空间离散化网格(单位：Å)

# 构造函数
```julia
Material(boxSizes::Vector{Float64}, inputGridVectors::Matrix{Float64}, atoms::Vector{Atom})
```

# 示例
```julia
# 创建材料系统
boxSizes = [100.0, 100.0, 50.0]  # 100×100×50 Å
gridVectors = [5.0 0.0 0.0; 0.0 5.0 0.0; 0.0 0.0 5.0]  # 5×5×5 Å 网格单元
material = Material(boxSizes, gridVectors, atoms)
```
"""
struct Material
    box::Box                          # 模拟盒子对象，定义材料的空间边界和周期性条件
    atoms::Vector{Atom}               # 原子向量，存储材料中的所有原子
    inputGridVectors::Matrix{Float64} # 3x3矩阵，输入网格向量，用于定义空间离散化网格
end

# 定义晶格点(LatticePoint)结构体，表示晶体结构中的固定位置
mutable struct LatticePoint
    index::Int64                      # 晶格点的唯一标识符
    type::Int64                       # 晶格点对应的原子类型，初始化后保持不变
    coordinate::Vector{Float64}       # 3维向量，晶格点的空间位置坐标(单位：Å)
    cellIndex::Tuple{Int64, Int64, Int64}  # 三元组，晶格点所在网格单元的索引(x,y,z)
    environment::Vector{Int64}        # 整数向量，存储邻近晶格点的索引，用于环境识别

    atomIndex::Int64                  # 占据该晶格点的原子索引，-1表示空位(vacancy)
end

# 定义邻近单元信息(NeighborCellInfo)结构体，存储相邻网格单元的信息
struct NeighborCellInfo
    index::NTuple{3, Int64}           # 三元组，邻近网格单元的索引(x,y,z)
    cross::NTuple{3, Int8}            # 三元组，周期性边界穿越标志：0-无穿越，1-正向穿越，-1-负向穿越
                                      # 例如：(0,0,1)表示在z方向正向穿越边界
end

"""
    Cell

表示空间离散化后的计算单元，用于空间分区和邻居搜索优化。

仅适用于正交盒子系统。

# 字段

## 基本字段
- `index::Tuple{Int64, Int64, Int64}`: 三元组，网格单元在网格中的索引(x,y,z)
- `atoms::Vector{Atom}`: 原子向量，存储位于该单元内的所有活动原子
- `latticePoints::Vector{LatticePoint}`: 晶格点向量，存储位于该单元内的所有晶格点
- `ranges::Matrix{Float64}`: 3x2矩阵，存储单元的空间边界范围[min, max]×3(单位：Å)
- `neighborCellsInfo::Array{NeighborCellInfo, 3}`: 3x3x3数组，存储所有邻近单元的信息（包括周期性边界）
- `isExplored::Bool`: 布尔标志，指示该单元在当前碰撞搜索中是否已被探索
- `atomicDensity::Float64`: 原子数密度，该单元内的平均原子密度(单位：原子/Å³)

## 动态加载系统字段
- `latticeAtoms::Vector{Atom}`: 原子向量，存储该单元内新加载的晶格原子（动态加载模式）
- `isLoaded::Bool`: 布尔标志，指示该单元的晶格原子是否已加载到内存
- `vacancies::Vector{Atom}`: 原子向量，存储该单元内的空位缺陷（静态和动态加载都使用）
- `isSavedLatticeRange::Bool`: 布尔标志，指示晶格范围是否已计算并保存
- `latticeRanges::Matrix{Int64}`: 3x2矩阵，存储该单元对应的晶格索引范围
- `isPushedNeighbor::Bool`: 布尔标志，指示邻近单元信息是否已初始化

# 构造函数
```julia
Cell(index, atoms, latticePoints, ranges, neighborCellsInfo, isExplored, atomicDensity)
```

# 示例
```julia
# 创建网格单元
index = (1, 1, 1)
atoms = Vector{Atom}()
latticePoints = Vector{LatticePoint}()
ranges = [0.0 5.0; 0.0 5.0; 0.0 5.0]  # 5×5×5 Å 范围
neighborInfo = Array{NeighborCellInfo, 3}(undef, 3, 3, 3)
cell = Cell(index, atoms, latticePoints, ranges, neighborInfo, false, 0.0)
```
"""
mutable struct Cell
    # 仅适用于正交盒子系统
    index::Tuple{Int64, Int64, Int64}  # 三元组，网格单元在网格中的索引(x,y,z)
    atoms::Vector{Atom}               # 原子向量，存储位于该单元内的所有原子
    latticePoints::Vector{LatticePoint}  # 晶格点向量，存储位于该单元内的所有晶格点
    ranges::Matrix{Float64}           # 3x2矩阵，存储单元的空间边界范围[min, max]×3(单位：Å)
    neighborCellsInfo::Array{NeighborCellInfo, 3}  # 3x3x3数组，存储所有邻近单元的信息
    isExplored::Bool                  # 布尔标志，指示该单元在当前碰撞搜索中是否已被探索
    atomicDensity::Float64            # 原子数密度，该单元内的平均原子密度(单位：原子/Å³)
    
    # 以下字段用于动态加载系统
    latticeAtoms::Vector{Atom}        # 原子向量，存储该单元内新加载的晶格原子
    isLoaded::Bool                    # 布尔标志，指示该单元的晶格原子是否已加载
    vacancies::Vector{Atom}           # 原子向量，存储该单元内的空位，也用于静态加载
    isSavedLatticeRange::Bool         # 布尔标志，指示晶格范围是否已计算并保存
    latticeRanges::Matrix{Int64}      # 3x2矩阵，存储该单元对应的晶格索引范围
    isPushedNeighbor::Bool            # 布尔标志，指示邻近单元信息是否已初始化
end

# 网格单元的构造函数，提供默认值初始化动态加载相关字段
function Cell(
    index::Tuple{Int64, Int64, Int64},
    atoms::Vector{Atom},
    latticePoints::Vector{LatticePoint},
    ranges::Matrix{Float64},
    #neighborCellsInfo::Dict{Vector{Int8}, NeighborCellInfo},  # 注释：旧版本的字典存储方式
    neighborCellsInfo::Array{NeighborCellInfo, 3},
    isExplored::Bool,
    atomicDensity::Float64)           
    latticeAtoms = Vector{Atom}()     # 初始化空的晶格原子向量
    isLoaded = false                  # 默认未加载晶格原子
    vacancies = Vector{Atom}()        # 初始化空的空位向量
    isSavedLatticeRange = false       # 默认未保存晶格范围
    latticeRanges = Matrix{Int64}(undef, 3, 2)  # 分配3x2矩阵用于存储晶格范围
    isPushedNeighbor = false          # 默认邻近信息未初始化
    return Cell(index, atoms, latticePoints, ranges, neighborCellsInfo, isExplored, atomicDensity, 
                    latticeAtoms, isLoaded, vacancies, isSavedLatticeRange, latticeRanges, isPushedNeighbor)     
end

# 宏定义，根据动态加载标志决定网格单元的存储类型
# 注意：此宏需要在有 Parameters 对象的环境中调用，或使用全局配置
# 为了兼容性，我们提供一个函数版本，接受 is_dynamic_load 参数
function get_cell_storage_type(is_dynamic_load::Bool)
    if is_dynamic_load
        return Dict{Tuple{Int64, Int64, Int64}, Cell}  # 动态加载使用字典存储，支持稀疏网格
    else
        return Array{Cell, 3}  # 静态加载使用三维数组存储，支持密集网格
    end
end

# 为了向后兼容，保留宏定义，但需要传入 parameters
macro cell_storage_type(parameters)
    return quote
        if $(parameters).is_dynamic_load
            Dict{Tuple{Int64, Int64, Int64}, Cell}
        else
            Array{Cell, 3}
        end
    end
end

"""
    Grid

管理整个模拟空间的离散化网格系统，用于空间分区和邻居搜索优化。

# 字段
- `cells::Union{Array{Cell, 3}, Dict{Tuple{Int64, Int64, Int64}, Cell}}`: 网格单元集合
  - 静态加载模式：使用 `Array{Cell, 3}` 三维数组存储所有单元
  - 动态加载模式：使用 `Dict{Tuple{Int64, Int64, Int64}, Cell}` 字典存储已加载的单元
- `vectors::Matrix{Float64}`: 3x3矩阵，单个网格单元的基向量(单位：Å)
- `sizes::Vector{Int64}`: 3维向量，网格在各维度的大小（单元数量）
- `cellVolume::Float64`: 单个网格单元的体积(单位：Å³)

# 注意
- `cells` 字段使用 `Union` 类型以支持静态（Array）和动态（Dict）两种存储方式
- 存储类型由 `Parameters.is_dynamic_load` 决定
- 动态加载模式可以显著减少内存使用，适用于大规模模拟

# 示例
```julia
# 网格通常由 Simulator 初始化时自动创建
# 静态加载：cells 为 Array{Cell, 3}
# 动态加载：cells 为 Dict{Tuple{Int64, Int64, Int64}, Cell}
```
"""
mutable struct Grid
    cells::Union{Array{Cell, 3}, Dict{Tuple{Int64, Int64, Int64}, Cell}}  # 网格单元集合，根据 is_dynamic_load 选择存储类型
    vectors::Matrix{Float64}          # 3x3矩阵，单个网格单元的基向量(单位：Å)
    sizes::Vector{Int64}              # 3维向量，网格在各维度的大小（单元数量）
    cellVolume::Float64               # 单个网格单元的体积(单位：Å³)
end 

# 定义类型相关常数(ConstantsByType)结构体，存储与原子类型相关的物理常数
struct ConstantsByType
    V_upterm::Dict{Vector{Int64}, Float64}  # 字典，键为[type_p, type_t]，值为势函数的上项常数
    a_U::Dict{Vector{Int64}, Float64}       # 字典，通用屏蔽长度常数(单位：Å)
    E_m::Dict{Int64, Float64}               # 字典，键为原子类型，值为电子阻止本领的特征能量(单位：eV)
    S_e_upTerm::Dict{Vector{Int64}, Float64}  # 字典，电子阻止本领的上项常数
    S_e_downTerm::Dict{Vector{Int64}, Float64}  # 字典，电子阻止本领的下项常数
    x_nl::Dict{Vector{Int64}, Float64}      # 字典，非局域能量损失的距离参数
    a::Dict{Vector{Int64}, Float64}         # 字典，指数衰减长度参数(单位：Å)
    Q_nl::Dict{Vector{Int64}, Float64}      # 字典，非局域能量损失系数
    Q_loc::Dict{Vector{Int64}, Float64}     # 字典，局域能量损失系数
    qMax::Dict{Vector{Int64}, Float64}      # 字典，最大碰撞参数限制(单位：Å)
    sigma::Dict{Int64, Float64}             # 字典，键为原子类型，值为热振动的均方根位移(单位：Å)
end

# 定义元素(Element)结构体，存储化学元素的物理性质
struct Element
    name::String                      # 元素名称（化学符号）
    radius::Float64                   # 原子半径(单位：Å)
    mass::Float64                     # 原子质量(单位：原子质量单位u)
    Z::Float64                        # 原子序数（核电荷数）
    dte::Float64                      # 位移阈值能量(单位：eV)
    bde::Float64                      # 结合能(单位：eV)
    alpha::Float64                    # 电子阻止本领的α参数（无量纲）
    beta::Float64                     # 电子阻止本领的β参数（无量纲）
end

"""
    Parameters

存储模拟的各种控制参数和设置。

这是 DISPLATH 的核心配置结构体，包含晶体结构定义、碰撞参数、能量参数、温度参数等所有模拟设置。

# 主要字段
- `primaryVectors::Matrix{Float64}`: 3×3矩阵，原胞基向量(单位：Å)
- `latticeRanges::Matrix{Int64}`: 3×2矩阵，晶格索引范围[min,max]×3
- `basisTypes::Vector{Int64}`: 原胞内各原子的类型标识符
- `basis::Matrix{Float64}`: N×3矩阵，原胞内原子的分数坐标
- `typeDict::Dict{Int64, Element}`: 原子类型定义字典
- `pMax::Float64`: 最大碰撞参数(单位：Å)
- `stopEnergy::Float64`: 停止能量阈值(单位：eV)
- `temperature::Float64`: 模拟温度(单位：K)
- `is_dynamic_load::Bool`: 是否启用动态加载模式

# 构造函数
```julia
Parameters(primaryVectors, latticeRanges, basisTypes, basis, pMax, 
           vacancyRecoverDistance, typeDict; kwargs...)
```

# 示例
```julia
# 创建参数对象
primaryVectors = [5.431 0.0 0.0; 0.0 5.431 0.0; 0.0 0.0 5.431]
latticeRanges = [0 10; 0 10; 0 10]
basisTypes = [1, 1]
basis = [0.0 0.0 0.0; 0.5 0.5 0.5]
typeDict = Dict(1 => Element("Si", 20.0, 10.0))
parameters = Parameters(primaryVectors, latticeRanges, basisTypes, basis,
                        1.8, 1.3, typeDict; stopEnergy=0.1, temperature=300.0)
```
"""
mutable struct Parameters
    # 必需参数 - 晶体结构定义
    primaryVectors::Matrix{Float64}   # 3x3矩阵，原胞基向量(单位：Å)
    primaryVectors_INV::Matrix{Float64}  # 3x3矩阵，原胞基向量的逆矩阵，自动计算非输入参数
    latticeRanges::Matrix{Int64}      # 3x2矩阵，晶格索引范围[min,max]×3
    basisTypes::Vector{Int64}         # 整数向量，原胞内各原子的类型标识符
    basis::Matrix{Float64}            # Nx3矩阵，原胞内原子的分数坐标
    
    # 文件路径参数
    θτRepository::String              # 字符串，θ和τ数据文件的存储目录路径
    
    # 碰撞参数
    pMax::Float64                     # 最大碰撞参数，碰撞搜索的截断距离(单位：Å)
    pMax_squared::Float64             # pMax的平方值，自动计算用于距离比较优化
    
    # 缺陷参数
    vacancyRecoverDistance_squared::Float64  # 空位恢复距离的平方值，自动计算用于距离比较
    
    # 元素类型定义
    typeDict::Dict{Int64, Element}    # 字典，键为类型标识符，值为Element结构体
    
    # 可选参数 - 边界条件和几何设置
    periodic::Vector{Bool}            # 3维布尔向量，各维度的周期性边界条件标志
    isOrthogonal::Bool                # 布尔标志，指示模拟盒子是否为正交
    isPrimaryVectorOrthogonal::Bool   # 布尔标志，指示原胞基向量是否正交，自动计算非输入参数
    
    # 能量和碰撞参数插值范围
    EPowerRange::StepRangeLen{Float64, Base.TwicePrecision{Float64}, Base.TwicePrecision{Float64}, Int64}  # 能量对数范围的步进序列(单位：log₁₀(eV))
    pPowerRange::StepRangeLen{Float64, Base.TwicePrecision{Float64}, Base.TwicePrecision{Float64}, Int64}  # 碰撞参数对数范围的步进序列(单位：log₁₀(Å))
    
    # 停止准则
    stopEnergy::Float64               # 停止能量，原子动能低于此值时停止运动(单位：eV)
    
    # 能量损失模型选项
    isNonQnl::Bool                    # 布尔标志，为true时忽略非局域能量损失
    
    # 温度相关参数
    DebyeTemperature::Float64         # 德拜温度(单位：K)，用于热振动计算
    isDumpInCascade::Bool             # 布尔标志，为true时在碰撞级联过程中输出原子位置
    DTEMode::Int64                    # 整数，位移阈值能量计算模式：1-直接，2-环境相关，4-自定义
    
    # SOAP描述符参数（注释掉，可能用于未来扩展）
    #soapParameters::Vector{Float64}
    
    # DTE数据文件
    DTEFile::String                   # 字符串，位移阈值能量数据文件路径
    
    # KMC模拟参数
    isKMC::Bool                       # 布尔标志，为true时启用动力学蒙特卡洛模拟
    nu_0_dict::Dict{Int64, Float64}   # 字典，键为原子类型，值为尝试频率(单位：Hz)
    temperature::Float64              # 模拟温度(单位：K)
    temperature_kb::Float64           # 以eV为单位的温度值(kB×T)，自动计算
    
    # 完美晶格环境索引
    perfectEnvIndex::Int64            # 整数，完美晶格环境的索引值
    
    # 辐照参数
    irrdiationFrequency::Float64      # 辐照频率(单位：Hz)，用于KMC模拟中的辐照事件
    
    # 动态加载参数
    is_dynamic_load::Bool             # 布尔标志，为true时启用动态加载模式（按需加载原子）
    nCascadeEveryLoad::Int64          # 整数，每进行多少次碰撞级联后清理晶格原子
    maxRSS::Int64                     # 整数，最大内存使用量(单位：kB)，用于内存管理
    
    # 非晶材料参数
    isAmorphous::Bool                 # 布尔标志，为true时模拟非晶材料
    amorphousHeight::Float64          # 浮点数，非晶区域的高度位置(单位：Å)
    
    # 随机数生成器参数
    random_seed::Int64                # 随机数生成器种子，用于可重现的随机数序列
    
    # 调试模式
    debugMode::Bool                   # 布尔标志，为true时启用调试输出和额外检查
end

# Parameters结构体的构造函数，用于创建和初始化模拟参数对象
function Parameters(
    # 必需参数 - 晶体结构和碰撞参数的基本定义
    primaryVectors::Matrix{Float64},           # 3x3矩阵，原胞基向量(单位：Å)，定义晶体周期性结构
    latticeRanges::Matrix{Int64},              # 3x2矩阵，晶格索引范围[min,max]×3，定义模拟晶格大小
    basisTypes::Vector{Int64},                 # 整数向量，原胞内各原子位置的类型标识符
    basis::Matrix{Float64},                    # Nx3矩阵，原胞内原子的分数坐标，定义晶体结构
    pMax::Float64,                             # 最大碰撞参数(单位：Å)，碰撞搜索的截断距离
    vacancyRecoverDistance::Float64,           # 空位恢复距离(单位：Å)，判定空位复合的距离阈值
    typeDict::Dict{Int64, Element};            # 字典，原子类型定义，键为类型ID，值为元素性质
    
    # 可选参数 - 模拟设置和高级选项，提供合理的默认值
    periodic::Vector{Bool} = [true, true, false],  # 3维布尔向量，各维度周期性边界条件，默认xy周期z开放
    isOrthogonal::Bool = true,                     # 布尔标志，模拟盒子是否正交，默认true简化计算
    θτRepository::String = ENV["ARCS_REPO"] * "/thetatau_repository/",  # θτ数据文件目录路径
    EPowerRange::StepRangeLen{Float64, Base.TwicePrecision{Float64}, Base.TwicePrecision{Float64}, Int64} = -1.0:0.045:8.0,  # 能量对数插值范围：10⁻¹到10⁸ eV，步长0.045
    pPowerRange::StepRangeLen{Float64, Base.TwicePrecision{Float64}, Base.TwicePrecision{Float64}, Int64} = -10.0:0.01:1.0,  # 碰撞参数对数插值范围：10⁻¹⁰到10¹ Å，步长0.01
    stopEnergy::Float64 = 0.1,                    # 停止能量阈值(单位：eV)，原子动能低于此值停止运动
    isNonQnl::Bool = false,                       # 布尔标志，为true时忽略非局域能量损失，仅静态加载有效
    DebyeTemperature::Float64 = 519.0,            # 德拜温度(单位：K)，用于热振动计算，默认石墨烯值
    isDumpInCascade::Bool = false,                # 布尔标志，为true时在碰撞级联过程中输出原子快照
    DTEMode::Int64 = 1,                           # 位移阈值能量计算模式：1-直接，2-环境相关，4-自定义
    #soapParameters::Vector{Float64} = [2.6, 8.0, 6.0],  # 注释：SOAP描述符参数[截断半径, n_max, l_max]
    DTEFile::String="",                           # 字符串，位移阈值能量数据文件路径，DTEMode=2时必需
    isKMC::Bool = false,                          # 布尔标志，为true时启用动力学蒙特卡洛模拟
    nu_0_dict::Dict{Int64, Float64} = Dict{Int64, Float64}(),  # 字典，KMC尝试频率(单位：Hz)，键为原子类型
    temperature::Float64 = 0.0,                   # 模拟温度(单位：K)，0K表示不考虑热振动
    perfectEnvIndex::Int64 = 0,                   # 整数，完美晶格环境的索引值，用于KMC环境识别
    irrdiationFrequency::Float64 = 0.0,           # 辐照频率(单位：Hz)，KMC模拟中辐照事件的发生频率
    is_dynamic_load::Bool = false,                 # 布尔标志，为true时启用动态加载模式（按需加载原子，节省内存）
    nCascadeEveryLoad = 100,                      # 整数，动态加载中每N次碰撞级联后执行内存清理
    maxRSS::Int = 20,                             # 整数，最大内存使用量(单位：GB)，用于动态内存管理
    isAmorphous = false,                          # 布尔标志，为true时模拟非晶材料结构
    amorphousLength::Float64 = -100.0,            # 浮点数，非晶区域长度(单位：Å)，从顶部开始计算
    random_seed::Int64 = 42,                      # 随机数生成器种子，用于可重现的随机数序列
    debugMode::Bool = false)                      # 布尔标志，为true时启用调试模式和额外检查
    
    # ========== 输入验证 ==========
    # 注意：异常类型在 logging.jl 中定义，通过 DISPLATH 模块导出
    # 由于 types.jl 在 DISPLATH 模块内部被 include，可以直接使用异常类型
    
    # 验证原胞基向量
    if size(primaryVectors) != (3, 3)
        throw(InvalidParameterError("primaryVectors", size(primaryVectors), "Must be a 3×3 matrix"))
    end
    det_primary = det(primaryVectors)
    if det_primary <= 0.0
        throw(InvalidParameterError("primaryVectors", det_primary, "Determinant must be positive (volume > 0)"))
    end
    
    # 验证晶格范围
    if size(latticeRanges) != (3, 2)
        throw(InvalidParameterError("latticeRanges", size(latticeRanges), "Must be a 3×2 matrix [min max]×3"))
    end
    for d in 1:3
        if latticeRanges[d, 1] >= latticeRanges[d, 2]
            throw(InvalidParameterError("latticeRanges[$d,:]", latticeRanges[d,:], "min must be < max"))
        end
    end
    
    # 验证基原子
    if length(basisTypes) != size(basis, 1)
        throw(InvalidParameterError("basisTypes", length(basisTypes), "Length must match number of basis atoms"))
    end
    if size(basis, 2) != 3
        throw(InvalidParameterError("basis", size(basis), "Must be N×3 matrix"))
    end
    
    # 验证碰撞参数
    if pMax <= 0.0
        throw(InvalidParameterError("pMax", pMax, "Must be positive"))
    end
    if vacancyRecoverDistance < 0.0
        throw(InvalidParameterError("vacancyRecoverDistance", vacancyRecoverDistance, "Must be non-negative"))
    end
    
    # 验证能量参数
    if stopEnergy <= 0.0
        throw(InvalidParameterError("stopEnergy", stopEnergy, "Must be positive"))
    end
    if temperature < 0.0
        throw(InvalidParameterError("temperature", temperature, "Must be non-negative"))
    end
    if DebyeTemperature < 0.0
        throw(InvalidParameterError("DebyeTemperature", DebyeTemperature, "Must be non-negative"))
    end
    
    # 验证类型字典
    if isempty(typeDict)
        throw(InvalidParameterError("typeDict", typeDict, "Cannot be empty"))
    end
    for (type_id, element) in typeDict
        if type_id <= 0
            throw(InvalidParameterError("typeDict keys", type_id, "Type IDs must be positive integers"))
        end
    end
    
    # 验证动态加载参数
    if nCascadeEveryLoad <= 0
        throw(InvalidParameterError("nCascadeEveryLoad", nCascadeEveryLoad, "Must be positive"))
    end
    if maxRSS <= 0
        throw(InvalidParameterError("maxRSS", maxRSS, "Must be positive"))
    end
    
    # 验证 DTE 模式
    if DTEMode == 2 && isempty(DTEFile)
        throw(InvalidParameterError("DTEFile", DTEFile, "Required when DTEMode == 2"))
    end
    
    # 验证必要的目录和文件存在性
    if !isdir(θτRepository)
        throw(InvalidParameterError("θτRepository", θτRepository, "Directory does not exist"))
    end
    
    # 参数预处理和计算 - 自动计算派生参数
    pMax_squared = pMax * pMax                    # 计算pMax的平方值，用于距离比较优化计算效率
    temperature_kb = temperature * 8.61733362E-5  # 将温度从K转换为eV：kB = 8.61733362×10⁻⁵ eV/K
    primaryVectors_INV = inv(primaryVectors)      # 计算原胞基向量的逆矩阵，用于分数坐标转换
    
    # 几何属性计算 - 判断原胞基向量是否正交以优化计算
    isPrimaryVectorOrthogonal = (primaryVectors[1,2] == 0.0 && primaryVectors[1,3] == 0.0 && 
                    primaryVectors[2,1] == 0.0 && primaryVectors[2,3] == 0.0 && 
                    primaryVectors[3,1] == 0.0 && primaryVectors[3,2] == 0.0)  # 检查非对角元素是否为零
    
    vacancyRecoverDistance_squared = vacancyRecoverDistance * vacancyRecoverDistance  # 计算距离平方值用于高效比较
    
    maxRSS *= 1048576  # 将内存限制从GB转换为kB：1 GB = 1024×1024 kB = 1,048,576 kB
    
    # 计算非晶区域高度：晶格顶部坐标减去非晶长度
    amorphousHeight = latticeRanges[3,2] * primaryVectors[3,3] - amorphousLength
    
    # 创建并返回Parameters结构体实例
    return Parameters(primaryVectors, primaryVectors_INV, latticeRanges, basisTypes, basis,
                      θτRepository, pMax, pMax_squared, vacancyRecoverDistance_squared, typeDict,
                      periodic, isOrthogonal, isPrimaryVectorOrthogonal,
                      EPowerRange, pPowerRange, stopEnergy, isNonQnl, DebyeTemperature, isDumpInCascade, 
                      DTEMode, 
                      #soapParameters,  # 注释掉的SOAP参数
                      DTEFile,
                      isKMC, nu_0_dict, temperature, temperature_kb, perfectEnvIndex, irrdiationFrequency,
                      is_dynamic_load, nCascadeEveryLoad, maxRSS, isAmorphous, amorphousHeight, 
                      random_seed, debugMode)
end 

"""
    CollisionParamsBuffers

用于临时存储碰撞计算中的中间变量，避免重复内存分配。

该结构体在碰撞级联模拟中用于批量处理多个目标原子的碰撞参数，提高计算效率。

# 字段
- `tanφList::Vector{Float64}`: 存储多个目标原子的tanφ值，φ为入射原子散射角
- `tanψList::Vector{Float64}`: 存储多个目标原子的tanψ值，ψ为目标原子散射角
- `E_tList::Vector{Float64}`: 存储多个目标原子的能量转移值(单位：eV)
- `x_pList::Vector{Float64}`: 存储多个目标原子的入射原子位移参数(单位：Å)
- `x_tList::Vector{Float64}`: 存储多个目标原子的目标原子位移参数(单位：Å)
- `Q_locList::Vector{Float64}`: 存储多个目标原子的局域能量损失(单位：eV)

# 构造函数
```julia
CollisionParamsBuffers()
```

# 使用
```julia
buffers = CollisionParamsBuffers()
EnsureCollisionCapacity!(buffers, 10)  # 确保有足够容量存储10个目标的参数
```
"""
mutable struct CollisionParamsBuffers
    tanφList::Vector{Float64}        # 存储多个目标原子的tanφ值，φ为入射原子散射角
    tanψList::Vector{Float64}        # 存储多个目标原子的tanψ值，ψ为目标原子散射角  
    E_tList::Vector{Float64}         # 存储多个目标原子的能量转移值(单位：eV)
    x_pList::Vector{Float64}         # 存储多个目标原子的入射原子位移参数(单位：Å)
    x_tList::Vector{Float64}         # 存储多个目标原子的目标原子位移参数(单位：Å)
    Q_locList::Vector{Float64}       # 存储多个目标原子的局域能量损失(单位：eV)
    
    # 构造函数，初始化空的缓冲区向量
    function CollisionParamsBuffers()
        return new(
            Vector{Float64}(), Vector{Float64}(), Vector{Float64}(),  # 初始化三个空向量
            Vector{Float64}(), Vector{Float64}(), Vector{Float64}()   # 初始化另外三个空向量
        )
    end
end

# 工作缓冲区结构体，管理模拟过程中的临时数据存储，支持多线程并行计算
mutable struct WorkBuffers
    coordinates::Vector{Vector{Float64}}        # 坐标缓冲区向量，每个线程一个3维坐标向量
    candidateTargets::Vector{Atom}              # 候选目标原子向量，存储当前搜索到的可能碰撞目标
    collisionParames::CollisionParamsBuffers    # 碰撞参数缓冲区，用于碰撞动力学计算
    threadCandidates::Vector{Vector{Atom}}      # 线程专用候选原子向量，每个线程一个，避免数据竞争
    threadRNG::Vector{StableRNGs.StableRNG}    # 线程本地随机数生成器，每个线程一个，确保线程安全
    
    # 构造函数，根据最大线程数和随机种子初始化缓冲区
    function WorkBuffers(max_threads::Int64=Threads.nthreads(), seed::Int64=42)
        # 为每个线程预分配3维坐标向量，避免动态内存分配
        coordinates = [Vector{Float64}(undef, 3) for _ in 1:max_threads]
        
        # 初始化候选目标原子向量，预分配容量以提高性能
        candidateTargets = Vector{Atom}()
        sizehint!(candidateTargets, 100)  # 预分配100个原子的容量
        
        # 创建碰撞参数缓冲区实例
        collisionParams = CollisionParamsBuffers()
        
        # 为每个线程创建专用的候选原子向量，避免多线程数据竞争
        threadCandidates = [Vector{Atom}() for _ in 1:max_threads]
        for tc in threadCandidates
            sizehint!(tc, 50)  # 每个线程预分配50个原子的容量
        end
        
        # 为每个线程创建独立的稳定随机数生成器，确保线程安全
        # 注意：StableRNG 在 DISPLATH 模块中导入
        threadRNG = [StableRNGs.StableRNG(seed + t) for t in 1:max_threads]
        
        # 返回初始化的WorkBuffers实例
        return new(coordinates, candidateTargets, 
                  collisionParams, threadCandidates, threadRNG)
    end
end

# 确保碰撞参数缓冲区有足够的容量存储N个目标的参数
function EnsureCollisionCapacity!(buffers::CollisionParamsBuffers, n::Int)
    if length(buffers.tanφList) < n
        # 调整所有参数向量的长度以适应N个目标原子
        resize!(buffers.tanφList, n)    # 调整tanφ向量容量
        resize!(buffers.tanψList, n)    # 调整tanψ向量容量  
        resize!(buffers.E_tList, n)     # 调整能量转移向量容量
        resize!(buffers.x_pList, n)     # 调整入射原子位移向量容量
        resize!(buffers.x_tList, n)     # 调整目标原子位移向量容量
        resize!(buffers.Q_locList, n)   # 调整局域能量损失向量容量
    end
end

# 清空工作缓冲区中的所有临时数据，为新的计算周期做准备
function ClearBuffers!(buffers::WorkBuffers)
    empty!(buffers.coordinate)          # 清空坐标缓冲区（注：字段名应为coordinates）
    empty!(buffers.targets)             # 清空目标原子缓冲区（注：字段名可能不匹配）
    empty!(buffers.candidateTargets)    # 清空候选目标原子向量
    
    # 清空所有线程专用的候选原子向量
    for tc in buffers.threadCandidates
        empty!(tc)  # 清空单个线程的候选原子向量
    end
end

# 主模拟器结构体，管理整个碰撞级联模拟的状态和数据
mutable struct Simulator
    atoms::Vector{Atom}                    # 原子向量，存储模拟系统中所有活动的原子
    latticePoints::Vector{LatticePoint}    # 晶格点向量，存储晶体结构的所有晶格位置
    box::Box                               # 模拟盒子，定义系统的空间边界和周期性
    grid::Grid                             # 空间网格，用于空间分区和邻居搜索优化
    maxAtomID::Int64                       # 当前最大原子ID，用于分配新原子的唯一标识符
    numberOfAtoms::Int64                   # 系统中活动原子的总数
    constantsByType::ConstantsByType       # 类型相关常数，存储原子类型相关的物理参数
    isStore::Bool                          # 存储标志，为true时系统状态已被保存用于恢复
    displacedAtoms::Vector{Int64}          # 位移原子索引向量，存储从晶格位置移出的原子ID
    numberOfAtomsWhenStored::Int64         # 存储时的原子数量，用于状态恢复的基准
    nCascade::Int64                        # 碰撞级联计数器，记录已完成的碰撞级联数量
    nCollisionEvent::Int64                 # 碰撞事件计数器，记录当前级联中的碰撞事件数
    exploredCells::Vector{Cell}            # 已探索单元向量，存储当前搜索中访问过的网格单元
    θFunctions::Dict{Vector{Int64}, Function}  # θ函数字典，键为[type_p,type_t]，值为散射角插值函数
    τFunctions::Dict{Vector{Int64}, Function}  # τ函数字典，键为[type_p,type_t]，值为时间参数插值函数
    uniformDensity::Float64                # 均匀密度，系统的平均原子数密度(单位：原子/Å³)
    #soap::PyObject                        # 注释：SOAP描述符计算对象，用于机器学习势
    environmentCut::Float64                # 环境截断距离(单位：Å)，用于定义局部原子环境
    DTEData::Vector{Vector{Float64}}       # DTE数据，存储不同环境下的位移阈值能量
    
    # 动力学蒙特卡洛(KMC)相关字段
    time::Float64                          # 模拟时间(单位：秒)，记录KMC模拟的累计时间
    frequency::Float64                     # 总频率(单位：Hz)，所有KMC事件频率之和
    frequencies::Vector{Float64}           # 频率向量，存储每个KMC事件的个体频率
    mobileAtoms::Vector{Atom}              # 可移动原子向量，存储当前可发生迁移的原子
    
    # 动态加载系统相关字段
    vacancies::Vector{Atom}                # 空位向量，存储系统中所有的空位缺陷
    numberOfVacancies::Int64               # 空位数量，系统中当前空位的总数
    maxVacancyID::Int64                    # 最大空位ID，用于分配新空位的唯一标识符
    minLatticeAtomID::Int64                # 最小晶格原子ID，动态加载中晶格原子的ID范围
    
    # 调试和开发相关字段
    debugAtoms::Vector{Atom}               # 调试原子向量，存储用于调试的额外原子信息
    parameters::Parameters                 # 模拟参数，存储所有的模拟设置和控制参数
    workBuffers::WorkBuffers               # 工作缓冲区，管理临时数据和计算缓存
end

# Simulator结构体的构造函数，初始化模拟器的主要组件和状态
function Simulator(box::Box, inputGridVectors::Matrix{Float64}, parameters::Parameters)
    # 初始化空间网格系统
    grid = CreateGrid(box, inputGridVectors, parameters)  # 根据盒子和输入网格向量创建空间分区网格
    
    # 初始化类型相关物理常数
    constantsByType = InitConstantsByType(parameters.typeDict, parameters)  # 计算原子类型相关的各种常数
    
    # 初始化散射角和时间参数插值函数
    θFunctions, τFunctions = InitθτFunctions(parameters, constantsByType)  # 加载或计算θ和τ的插值函数
    
    # 初始化SOAP描述符（注释状态）
    #soap = InitSoap(parameters)
    
    # 根据DTE模式加载位移阈值能量数据
    if parameters.DTEMode == 2
        environmentCut, DTEData = LoadDTEData(parameters)  # 模式2：从文件加载环境相关的DTE数据
    else
        environmentCut, DTEData = -1.0, Vector{Vector{Float64}}()  # 其他模式：使用默认值
    end
    
    # 初始化KMC相关字段
    time = 0.0                           # 初始模拟时间为零
    frequency = 0.0                      # 初始总频率为零
    frequencies = Vector{Float64}()      # 空频率向量
    mobileAtoms = Vector{Atom}()         # 空可移动原子向量
    
    # 初始化动态加载相关字段
    vacancies = Vector{Atom}()           # 空位向量
    nCollisionEvent = 0                  # 碰撞事件计数器归零
    numberOfVacancies = 0                # 初始空位数量为零
    maxVacancyID = 1E6                   # 设置空位ID的起始值（与原子ID区分）
    minLatticeAtomID = 0                 # 晶格原子ID起始值
    
    # 初始化调试相关字段
    debugAtoms = Atom[]                  # 空调试原子数组
    workBuffers = WorkBuffers(Threads.nthreads(), parameters.random_seed)  # 创建工作缓冲区实例，传入随机种子
    
    # 计算系统的均匀原子数密度：原胞内原子数除以原胞体积
    uniformDensity = length(parameters.basisTypes) / (parameters.primaryVectors[1,1] * parameters.primaryVectors[2,2] * parameters.primaryVectors[3,3])
    
    # 创建并返回完整的Simulator实例
    return Simulator(Vector{Atom}(), Vector{LatticePoint}(),  # 空的原子和晶格点向量
                     box, grid,                              # 几何结构组件
                     0, 0,                                   # 初始原子ID和数量为零
                     constantsByType,                        # 类型相关常数
                     false, Vector{Int64}(), 0,              # 存储状态相关字段初始值
                     0, nCollisionEvent,                     # 级联和碰撞计数器
                     Vector{Cell}(),                         # 空的已探索单元向量
                     θFunctions, τFunctions,                 # 散射插值函数
                     uniformDensity,                         # 均匀密度
                     #soap,                                 # 注释的SOAP对象
                     environmentCut, DTEData,               # DTE相关数据
                     time, frequency, frequencies, mobileAtoms,  # KMC相关字段
                     vacancies, numberOfVacancies, maxVacancyID, minLatticeAtomID,  # 动态加载字段
                     debugAtoms,                             # 调试字段
                     parameters,                             # 模拟参数
                     workBuffers)                            # 工作缓冲区
end