# 功能：统计模拟系统中的空位数量
function CountVacancies(simulator::Simulator)
    nVacancies = 0  # 计数器初始化：空位数量初始化为0
    for latticePoint in simulator.latticePoints  # 格点遍历：遍历所有晶格格点
        if latticePoint.atomIndex == -1  # 空位判断：atomIndex == -1表示该格点为空位，约定：-1表示没有原子占据该格点
            nVacancies += 1  # 计数增加：发现空位时计数器加1
        end
    end
    return nVacancies  # 返回结果：空位总数
end

# 功能：提取所有空位格点的列表
function ExtractVacancyLattices(simulator::Simulator)
    vacancyLattices = []  # 列表初始化：创建空列表存储空位格点，类型推断：Julia会自动推断为Vector{Any}
    for latticePoint in simulator.latticePoints  # 格点遍历：遍历所有晶格格点
        if latticePoint.atomIndex == -1  # 空位判断：检查格点是否为空位
            push!(vacancyLattices, latticePoint)  # 添加空位：将空位格点添加到列表中，函数调用：Julia内置的push!函数
        end
    end
    return vacancyLattices  # 返回结果：包含所有空位格点的数组
end