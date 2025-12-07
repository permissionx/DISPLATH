# 注释掉的Python库导入，可能用于未来扩展的原子环境分析
#using PyCall
#@pyimport ase  # 原子模拟环境(ASE)库
#@pyimport dscribe.descriptors as descriptors  # 结构描述符库

# 获取原子的位移阈值能量(Displacement Threshold Energy)
function GetDTE(atom::Atom, simulator::Simulator)
    if simulator.parameters.DTEMode == 1  # 模式1：直接使用原子自身的DTE值
        return atom.dte  # 直接返回原子结构体中存储的位移阈值能量(单位：eV)
    elseif simulator.parameters.DTEMode == 2  # 模式2：基于原子周围环境计算DTE
        return GetDTEByEnvironment(atom, simulator)  # 调用环境相关的DTE计算函数
    #elseif simulator.parameters.DTEMode == 3   # 模式3：使用SOAP描述符计算DTE（注释掉）
    #    return GetDTEBySoap(atom, simulator)  # 调用基于SOAP的DTE计算函数
    elseif simulator.parameters.DTEMode == 4  # 模式4：自定义DTE计算方法
        return GetDTECustom(atom, simulator)  # 调用自定义DTE计算函数
    end
end

# 获取原子的结合能(Binding Energy)
function GetBDE(atom::Atom, simulator::Simulator)  # BDE: binding energy
    if simulator.parameters.DTEMode == 1  # 模式1：直接使用原子自身的BDE值
        return atom.bde  # 直接返回原子结构体中存储的结合能(单位：eV)
    elseif simulator.parameters.DTEMode == 2  # 模式2：基于原子周围环境计算BDE
        return GetBDEByEnvironment(atom, simulator)  # 调用环境相关的BDE计算函数
    #elseif simulator.parameters.DTEMode == 3   # 模式3：使用SOAP描述符计算BDE（注释掉）
    #    return GetBDEBySoap(atom, simulator)  # 调用基于SOAP的BDE计算函数
    elseif simulator.parameters.DTEMode == 4  # 模式4：自定义BDE计算方法
        return GetBDECustom(atom, simulator)  # 调用自定义BDE计算函数
    end
end

# 获取晶格点的位移阈值能量
function GetDTE(latticePoint::LatticePoint, simulator::Simulator)
    index = GetEnvironmentIndex(latticePoint, simulator)  # 计算晶格点的环境索引
    # 根据晶格点类型和环境索引从DTE数据中查找对应的位移阈值能量
    return simulator.DTEData[latticePoint.type][index]
end

# 基于原子环境计算位移阈值能量
function GetDTEByEnvironment(atom::Atom, simulator::Simulator)
    # 检查原子是否在晶格点上且晶格点类型与原子类型一致
    if atom.latticePointIndex != -1 && simulator.latticePoints[atom.latticePointIndex].type == atom.type
        # 如果原子在正确类型的晶格点上，获取对应晶格点的DTE
        return GetDTE(simulator.latticePoints[atom.latticePointIndex], simulator)
    else
        return 0.1  # 默认值，当原子不在晶格点或类型不匹配时使用(单位：eV)
    end
end

# 基于原子环境计算结合能
function GetBDEByEnvironment(atom::Atom, simulator::Simulator)
    # 在当前实现中，环境相关的结合能等于位移阈值能量
    # 这可能是一个简化假设，或者后续会修改为独立的计算
    return GetDTE(atom, simulator)
end

# 注释掉的SOAP描述符初始化函数
#function InitSoap(parameters::Parameters)
#    soapParameters = parameters.soapParameters  # 获取SOAP参数
#    soap = descriptors.SOAP(  # 创建SOAP描述符对象
#        species = [type.name for type in values(parameters.typeDict)],  # 元素种类列表
#        periodic = true,  # 启用周期性边界条件
#        r_cut = soapParameters[1],  # 截断半径(单位：Å)
#        n_max = Int(soapParameters[2]),  # 径向基函数的最大数量
#        l_max = Int(soapParameters[3])   # 球谐函数的最大角动量量子数
#    )
#    return soap  # 返回初始化的SOAP描述符
#end

# 获取原子周围邻居原子的坐标和元素名称数组
function GetNeighborArray(atom::Atom, simulator::Simulator)
    coordinates = atom.coordinate'  # 将原子坐标转置为行向量格式
    elementNames = Vector{String}([simulator.parameters.typeDict[atom.type].name])  # 初始化元素名称向量，包含中心原子的元素
    
    typeDict = simulator.parameters.typeDict  # 获取类型字典引用
    grid = simulator.grid  # 获取网格引用
    theCell = GetCell(grid, atom.cellIndex)  # 获取原子所在的网格单元
    
    # 遍历所有邻近单元信息
    for neighborInfo in theCell.neighborCellsInfo
        index = neighborInfo.index  # 获取邻近单元索引
        theCell = GetCell(grid, index)  # 获取邻近单元对象
        
        # 遍历邻近单元中的所有原子
        for cellAtom in cell.atoms
            if cellAtom.index != atom.index  # 排除中心原子自身
                coordinate = cellAtom.coordinate'  # 获取邻居原子坐标并转置为行向量
                coordinates = vcat(coordinates, coordinate)  # 垂直拼接坐标到坐标数组
                elementNames = push!(elementNames, typeDict[cellAtom.type].name)  # 添加邻居原子元素名称
            end
        end
    end
    
    return coordinates, elementNames  # 返回邻居原子坐标数组和元素名称向量
end

# 注释掉的SOAP描述符创建函数
#function CreateSoap(atom::Atom, simulator::Simulator)
#    # 构建模拟盒子的晶格向量
#    cell = [simulator.box.vectors[i,i] for i in 1:3]
#    
#    # 获取原子周围邻居的坐标和元素名称
#    coordinates, elementNames = GetNeighborArray(atom, simulator)
#    
#    pbc = true  # 设置周期性边界条件为真
#    
#    # 创建ASE原子对象
#    atoms = ase.Atoms(elementNames, coordinates, cell = cell, pbc = pbc)
#    
#    soap = simulator.soap  # 获取SOAP描述符对象
#
#    # 创建SOAP描述符，以原子位置0为中心
#    descriptor = soap.create(atoms, centers=[0])
#    
#    return soap, descriptor  # 返回SOAP对象和描述符
#end