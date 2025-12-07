function Parameters(pMax::Float64, vacancyRecoverDistance::Float64; kwargs...)
    # non lattice info
    primaryVectors = [1.0 0.0 0.0; 0.0 1.0 0.0; 0.0 0.0 1.0]
    latticeRanges = [0 1; 0 1; 0 1]
    basis = [0.0 0.0 0.0]
    basisTypes = [1]  
    typeDict = Dict{Int64, Element}()
    parameters = Parameters(primaryVectors, latticeRanges, basisTypes, basis, pMax, vacancyRecoverDistance, typeDict; kwargs...)
    return parameters
end

function Material(
    primaryVectors::Matrix{Float64},
    latticeRanges::Matrix{Int64},
    basisTypes::Vector{Int64},
    basis::Matrix{Float64},
    typeDict::Dict{Int64, Element},
    boxSizes::Vector{Int64}, 
    inputGridVectors::Matrix{Float64},
    parameters::Parameters)
    # 输入验证
    if length(boxSizes) != 3
        throw(InvalidParameterError("boxSizes", length(boxSizes), "Must be a 3-element vector"))
    end
    if any(s -> s <= 0, boxSizes)
        throw(InvalidParameterError("boxSizes", boxSizes, "All sizes must be positive"))
    end
    if size(inputGridVectors) != (3, 3)
        throw(InvalidParameterError("inputGridVectors", size(inputGridVectors), "Must be a 3×3 matrix"))
    end
    
    atoms = _PassingLatticeParamtersAndCreateAtoms(primaryVectors, latticeRanges, basisTypes, basis, typeDict, parameters)
    box = CreateBoxByPrimaryVectors(parameters.primaryVectors, boxSizes)
    material = Material(box, atoms, inputGridVectors)
    return material
end

function Material(
    primaryVectors::Matrix{Float64},
    latticeRanges::Matrix{Int64},
    basisTypes::Vector{Int64},
    basis::Matrix{Float64},
    typeDict::Dict{Int64, Element},
    boxVectors::Matrix{Float64}, 
    inputGridVectors::Matrix{Float64},
    parameters::Parameters)
    # 输入验证
    if size(inputGridVectors) != (3, 3)
        throw(InvalidParameterError("inputGridVectors", size(inputGridVectors), "Must be a 3×3 matrix"))
    end
    
    atoms = _PassingLatticeParamtersAndCreateAtoms(primaryVectors, latticeRanges, basisTypes, basis, typeDict, parameters)
    box = Box(boxVectors)  # Box 构造函数内部已有验证
    material = Material(box, atoms, inputGridVectors)
    return material
end

function Material(fileName::String, typeDict::Dict{Int64, Element}, inputGridVectors::Matrix{Float64}, parameters::Parameters; replicate::Vector{Int64} = [1,1,1])
    # In this mode, primaryVectors, latticeRanges, basisTypes, basis in parameters are expired.
    if parameters.is_dynamic_load
        throw(SimulationError("Material construction", "Material from date file is not supported in dynamic load mode."))
    end 
    parameters.typeDict = typeDict
    box, atoms = LoadAtomsAndBoxFromDataFile(fileName; replicate=replicate)
    material = Material(box, atoms, inputGridVectors)
    return material
end

function Simulator(material::Material, parameters::Parameters)
    simulator = Simulator(material.box, material.atoms, material.inputGridVectors, parameters)
    return simulator
end


# utils
function _PassingLatticeParamtersAndCreateAtoms(
    primaryVectors::Matrix{Float64},
    latticeRanges::Matrix{Int64},
    basisTypes::Vector{Int64},
    basis::Matrix{Float64},
    typeDict::Dict{Int64, Element},
    parameters::Parameters)
    parameters.primaryVectors = primaryVectors
    parameters.latticeRanges = latticeRanges
    parameters.basisTypes = basisTypes
    parameters.basis = basis
    parameters.typeDict = typeDict
    if !parameters.is_dynamic_load
        atoms = CreateAtomsByPrimaryVectors(parameters)
    else
        atoms = Atom[]
    end
    return atoms
end

# =====================================================================
# 通用辐照函数
# =====================================================================

"""
    Irradiation(simulator::Simulator, energy::Float64, ion_type::Int64, 
                ion_position::Vector{Float64}, ion_direction::Vector{Float64};
                restore::Bool=true, return_defects::Bool=false)

执行单次离子辐照模拟。

# 参数
- `simulator::Simulator`: 模拟器对象
- `energy::Float64`: 入射离子能量（单位：eV）
- `ion_type::Int64`: 离子类型（对应 typeDict 中的键）
- `ion_position::Vector{Float64}`: 离子初始位置（3D坐标，单位：Å）
- `ion_direction::Vector{Float64}`: 离子初始速度方向（3D单位向量）

# 关键字参数
- `restore::Bool=true`: 是否在辐照前恢复模拟器状态（静态加载模式需要先 Save!）
- `return_defects::Bool=false`: 是否返回缺陷统计信息

# 返回值
- 如果 `return_defects=false`: 返回空位数量（Int64）
- 如果 `return_defects=true`: 返回 `(nV, vacancies, interstitials)` 元组

# 示例
```julia
# 简单使用
nV = Irradiation(simulator, 1000.0, 2, [0.0, 0.0, 50.0], [0.0, 0.0, -1.0])

# 获取详细缺陷信息
nV, Vs, Is = Irradiation(simulator, 1000.0, 2, [0.0, 0.0, 50.0], 
                         [0.0, 0.0, -1.0], return_defects=true)
```
"""
function Irradiation(simulator::Simulator, energy::Float64, ion_type::Int64,
                     ion_position::Vector{Float64}, ion_direction::Vector{Float64};
                     restore::Bool=true, return_defects::Bool=false)
    # 输入验证
    if length(ion_position) != 3
        throw(InvalidParameterError("ion_position", length(ion_position), "Must be a 3-element vector"))
    end
    if length(ion_direction) != 3
        throw(InvalidParameterError("ion_direction", length(ion_direction), "Must be a 3-element vector"))
    end
    if energy <= 0.0
        throw(InvalidParameterError("energy", energy, "Must be positive"))
    end
    if !haskey(simulator.parameters.typeDict, ion_type)
        throw(InvalidParameterError("ion_type", ion_type, "Not found in typeDict"))
    end
    
    # 恢复模拟器状态（如果需要）
    if restore
        Restore!(simulator)
    end
    
    # 创建入射离子
    ion = Atom(ion_type, ion_position, simulator.parameters)
    SetVelocityDirection!(ion, ion_direction)
    SetEnergy!(ion, energy)
    push!(simulator, ion)
    
    # 执行碰撞级联
    Cascade!(ion, simulator)
    
    # 统计缺陷
    if return_defects
        _, Vs, Is = DefectStatics(simulator)
        return length(Vs), Vs, Is
    else
        _, Vs = DefectStatics(simulator)
        return length(Vs)
    end
end

"""
    Irradiation(simulator::Simulator, energy::Float64, ion_type::Int64;
                position_generator::Function, direction_generator::Function,
                restore::Bool=true, return_defects::Bool=false)

使用生成器函数创建离子位置和方向的辐照函数。

# 参数
- `simulator::Simulator`: 模拟器对象
- `energy::Float64`: 入射离子能量（单位：eV）
- `ion_type::Int64`: 离子类型

# 关键字参数
- `position_generator::Function`: 生成离子位置的函数 `() -> Vector{Float64}`
- `direction_generator::Function`: 生成离子方向的函数 `() -> Vector{Float64}`
- `restore::Bool=true`: 是否在辐照前恢复模拟器状态
- `return_defects::Bool=false`: 是否返回缺陷统计信息

# 示例
```julia
# 使用随机位置生成器
nV = Irradiation(simulator, 1000.0, 2;
                 position_generator=() -> RandomInSquare(41.6, 48.19, simulator) + [0.1, 0.1, 33.0],
                 direction_generator=() -> [0.0, 0.0, -1.0])
```
"""
function Irradiation(simulator::Simulator, energy::Float64, ion_type::Int64;
                     position_generator::Function, direction_generator::Function,
                     restore::Bool=true, return_defects::Bool=false)
    ion_position = position_generator()
    ion_direction = direction_generator()
    return Irradiation(simulator, energy, ion_type, ion_position, ion_direction;
                      restore=restore, return_defects=return_defects)
end

