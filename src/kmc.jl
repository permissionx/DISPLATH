#KMC（Kinetic Monte Carlo）动力学蒙特卡洛模拟

function UpdateEvents!(latticePointIndexs::Set{Int64}, simulator::Simulator) #更新指定格点集合的事件信息,latticePointIndexs::Set{Int64}：格点索引集合，使用Set确保唯一性,simulator::Simulator：主模拟器对象
    if simulator.parameters.isKMC #条件检查：仅在启用KMC模式时执行
        for latticePointIndex in latticePointIndexs #循环：遍历所有需要更新的格点
            latticePoint = simulator.latticePoints[latticePointIndex] #数据访问：从模拟器的格点数组中获取当前格点对象
            if latticePoint.atomIndex != -1 #条件检查：确保格点上有原子（atomIndex = -1 表示空位）
                envIndex = GetEnvironmentIndex(latticePoint, simulator) #函数调用：获取格点的环境索引
                atom = simulator.atoms[latticePoint.atomIndex] #数据访问：通过格点的atomIndex获取对应的原子对象
                if envIndex == simulator.parameters.perfectEnvIndex #条件检查：判断是否为完美晶格环境
                    DeleteAtomEvents!(simulator, atom) #函数调用：删除完美环境中原子的事件
                else
                    if atom.eventIndex != -1 #条件检查：原子已有事件记录（eventIndex = -1 表示无事件）
                        simulator.frequency -= atom.frequency #频率更新：从总频率中减去原子的旧频率
                        SetAtomEvents!(atom, latticePoint, simulator) #函数调用：重新设置原子的事件参数
                        simulator.frequency += atom.frequency #频率更新：添加新频率到总频率
                        simulator.frequencies[atom.eventIndex] = atom.frequency #数组更新：更新频率数组中对应位置的值
                    else
                        SetAtomEvents!(atom, latticePoint, simulator) #函数调用：设置原子的事件参数
                        AppendAtomEvents!(simulator, atom) #函数调用：将原子添加到事件列表
                    end
                end
            end
        end
    end
end

function SetAtomEvents!(atom::Atom, latticePoint::LatticePoint, simulator::Simulator) #设置原子的迁移事件参数
    finalLatticePointEnvIndexs, finalLatticePointIndexs = GetFinalLatticePointInfo(latticePoint, simulator) #函数调用：获取可迁移的目标格点信息
    energyBarriers = [simulator.energyBarrierDict[atom.type][[envIndex, i]] for i in finalLatticePointEnvIndexs] #能量势垒计算：根据原子类型和环境索引获取能垒
    #ComputeKMCProbability函数未定义
    atom.frequencies = [ComputeKMCProbability(type, e, simulator.temperature) for e in energyBarriers] #跃迁频率计算：根据能垒和温度计算阿伦尼乌斯跃迁频率,
    atom.frequency = sum(atom.frequencies) #总频率计算：原子所有可能跃迁路径的频率之和
    atom.finalLatticePointIndexs = finalLatticePointIndexs #目标格点存储：保存可迁移的目标格点索引
end

function AppendAtomEvents!(simulator::Simulator, atom::Atom) #将原子添加到可迁移原子列表
    push!(simulator.frequencies, atom.frequency) #数组操作：将原子频率添加到总频率数组
    push!(simulator.mobileAtoms, atom) #数组操作：将原子添加到可迁移原子数组
    simulator.maxEventIndex += 1 #索引管理：更新最大事件索引
    atom.eventIndex = simulator.maxEventIndex  #索引分配：为原子分配事件索引
    simulator.frequency += atom.frequency #总频率更新：累加到模拟器的总跃迁频率
end 

function GetFinalLatticePointInfo(latticePoint::LatticePoint, simulator::Simulator) #获取当前格点周围可迁移的目标格点
    environment = latticePoint.environment #环境获取：当前格点的邻居环境
    latticePoints = simulator.latticePoints #格点数组：所有格点的引用
    finalLatticePointEnvIndexs = Vector{Int64}() #初始化：创建空数组存储目标格点环境索引
    finalLatticePointIndexs = Vector{Int64}() #初始化：创建空数组存储目标格点实际索引
    for i in 1:length(environment) #循环：遍历所有邻居格点
        neighborLatticePoint = latticePoints[environment[i]] #邻居获取：通过环境索引获取邻居格点
        if neighborLatticePoint.atomIndex == -1 &&  latticePoint.tpye == atom.type #条件检查：邻居为空位且类型匹配（注意：这里有拼写错误`tpye`应为`type`，且`atom`未定义）
            push!(finalLatticePointEnvIndexs, 2^(i-1))    #环境索引计算：使用二进制位表示环境状态，第i个位置为1表示该邻居位置为空位
            push!(finalLatticePointIndexs, environment[i]) #索引存储：保存目标格点的实际索引
        end
    end
    return finalLatticePointEnvIndexs, finalLatticePointIndexs #返回值：环境索引数组和实际索引数组
end

function SetTemperature!(simulator::Simulator, temperature::Float64) #设置模拟温度
    temperature_kb = temperature * 8.61733362E-5 #单位转换：将开尔文温度转换为eV单位（玻尔兹曼常数$k_B = 8.61733362\times10^{-5}$ eV/K）
    simulator.parameters.temperature_kb = temperature_kb #参数更新：存储转换后的温度值
end

function ComputeProbability(type::Int64, energyBarrier::Float64, simulator::Simulator) #计算跃迁概率
    return simulator.parameters.nu_0_dict[type] * exp(-energyBarrier / simulator.parameters.temperature_kb) #阿伦尼乌斯公式：$ν = ν_0 \exp\left(-\frac{E_a}{k_B T}\right)$，其中nu_0_dict[type]：指前因子，与原子类型相关，energyBarrier：能垒$E_a$，temperature_kb：温度$k_B T$
end

function GetRandomMobileAtomIndex(simulator::Simulator) #使用拒绝采样方法随机选择迁移原子
    randomNumber = rand() * (simulator.frequency + simulator.parameters.irrdiationFrequency) #随机数生成：在[0, 总频率+辐照频率]区间均匀采样
    cumulativeProbability = 0.0 #初始化：累积概率初始化为0
    for i in 1:length(simulator.frequencies) #循环：遍历所有频率
        cumulativeProbability += simulator.frequencies[i] #概率累积：累加当前频率
        if cumulativeProbability >= randomNumber #条件检查：累积概率超过随机数
            return i #返回：对应的原子索引
        end
    end
    return -1 #特殊情况：未找到对应原子，可能选择辐照事件
end

function GetRandomFinalLatticePointIndex(atom::Atom) #为选定原子随机选择目标迁移格点
    randomNumber = rand() * atom.frequency #随机数生成：在[0, 原子总频率]区间采样
    cumulativeProbability = 0.0 #初始化：累积概率初始化为0
    for i in 1:length(atom.frequencies) #循环：遍历原子的所有跃迁频率
        cumulativeProbability += atom.frequencies[i] #概率累积：累加当前跃迁频率
        if randomNumber <= cumulativeProbability #条件检查：随机数小于等于累积概率
            return atom.finalLatticePointIndexs[i] #返回：对应的目标格点索引
        end
    end
end

function ProcessAnEvent!(simulator::Simulator) #处理一个KMC事件
    idx = GetRandomMobileAtomIndex(simulator) #原子选择：随机选择要迁移的原子
    if idx > 0 #条件检查：选择了有效的原子索引
        atom = simulator.mobileAtoms[idx] #数据访问：获取选中的原子对象
        latticePoint = simulator.latticePoints[GetRandomFinalLatticePointIndex(atom)] #目标格点选择：随机选择目标迁移格点
        Migrate!(atom, latticePoint, simulator) #函数调用：执行原子迁移操作
    else
        Irradiate!(simulator) #辐照事件：如果没有原子迁移，则执行辐照
    end
    ElapseTime!(simulator) #时间推进：根据总频率推进模拟时间
end

function DeleteAtomEvents!(simulator::Simulator, atom::Atom) #从事件列表中删除原子
    index = atom.eventIndex #索引获取：原子的当前事件索引
    simulator.frequency -= simulator.frequencies[index] #总频率更新：减去该原子的贡献
    deleteat!(simulator.frequencies, index) #数组操作：从频率数组中删除对应元素
    deleteat!(simulator.mobileAtoms, index) #数组操作：从可迁移原子数组中删除对应元素
    atom.eventIndex = -1 #状态重置：标记原子无事件
end 

function ElapseTime!(simulator::Simulator) #推进KMC模拟时间
    time = -1/simulator.frequency * log(rand()) #时间步长计算：使用泊松过程的时间步长公式，$\Delta t = -\frac{\ln(\xi)}{\sum \nu_i}$，其中$\xi \sim U(0,1)$
    simulator.time += time #时间累积：更新模拟器的总时间
end

function Migrate!(atom::Atom, latticePoint::LatticePoint, simulator::Simulator) #执行原子迁移操作
    oldLatticePoint = simulator.latticePoints[atom.latticePointIndex] #原格点获取：原子当前所在的格点
    LeaveLatticePoint!(atom, simulator; isUpdateEnv = false) #函数调用：原子离开当前格点
    SetOnLatticePoint!(atom, latticePoint, simulator; isUpdateEnv = false) #函数调用：原子占据新格点
    latticePointIndexs = Set([oldLatticePoint.environment; latticePoint.environment; latticePoint.index]) #受影响格点：收集需要更新事件的所有格点（旧环境+新环境+新格点）
    UpdateEvents!(latticePointIndexs, simulator) #事件更新：递归调用更新受影响格点的事件
end