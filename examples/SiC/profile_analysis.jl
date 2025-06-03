using Profile
using BenchmarkTools
using InteractiveUtils

# 包含主程序
home = "/beegfs/home/xuke/Researches/Irradiation_Li-Tianzhao/4.DISPLATH/DISPLATH/"
include(home * "src/main.jl")

println("=== DISPLATH 性能分析脚本 ===\n")

# 初始化参数 (复制from run.jl)
a = 4.36
primaryVectors = [a 0.0 0.0; 0.0 a 0.0; 0.0 0.0 a]
boxSizes = [50, 50, 1400]
inputGridVectors = [a*2.1 0.0 0.0; 0.0 a*2.1 0.0; 0.0 0.0 a*2.1]
latticeRanges = [0 50; 0 50; 2 1200]
basis = [0.0 0.0 0.0; 0.25 0.25 0.25]
basisTypes = [1, 2]

θτRepository = home * "thetatau_repository/"
pMax = 4.0
vacancyRecoverDistance = 4.0
typeDict = Dict(
    1 => Element("Si", 43.0, 22.0),  
    2 => Element("C", 40.0, 20.0),
    3 => Element("N", 1.0, 1.0)    
)
temperature = 300.0
DebyeTemperature = 490.0
parameters = Parameters(θτRepository, pMax, vacancyRecoverDistance, typeDict; 
                        temperature=temperature, DebyeTemperature=DebyeTemperature)

println("1. 初始化模拟器...")
@time simulator = Simulator(primaryVectors, boxSizes, inputGridVectors, latticeRanges, basis, basisTypes, parameters)

println("\n2. 保存初始状态...")
@time Save!(simulator)

println("\n=== 详细性能分析 ===")

# 创建测试离子
energy = 10000.0
ionPosition = [110.0, 110.0, 5240.0] 
ion = Atom(3, ionPosition, parameters)
SetVelocityDirection!(ion, [0.,0.,-1.])
SetEnergy!(ion, energy)
push!(simulator, ion)

println("\n3. 单次Cascade性能分析:")
println("3a. @time 分析:")
@time Cascade!(ion, simulator)

# 重置状态
Restore!(simulator)

# 创建新离子用于Profile分析
ion2 = Atom(3, ionPosition, parameters)
SetVelocityDirection!(ion2, [0.,0.,-1.])
SetEnergy!(ion2, energy)
push!(simulator, ion2)

println("\n3b. Profile 分析:")
Profile.clear()
@profile Cascade!(ion2, simulator)

println("\n=== Profile 结果 ===")
Profile.print(maxdepth=15, mincount=10)

println("\n=== 保存Profile结果到文件 ===")
open("profile_results.txt", "w") do io
    Profile.print(io, maxdepth=20, mincount=5)
end
println("详细Profile结果已保存到 profile_results.txt")

# 重置状态用于微基准测试
Restore!(simulator)

println("\n4. 关键函数微基准测试:")

# 测试LoadCell函数
println("4a. LoadCell函数:")
cell = simulator.cellGrid.cells[1,1,1]
cell.isLoaded = false  # 重置加载状态
@benchmark LoadCell($cell, $simulator) setup=(cell.isLoaded = false)

# 测试TemperatureToSigma函数
println("\n4b. TemperatureToSigma函数:")
@benchmark TemperatureToSigma(300.0, 490.0, 28.0)

# 测试原子创建
println("\n4c. Atom创建:")
@benchmark Atom(1, [0.0, 0.0, 0.0], $parameters)

# 测试矩阵计算
println("\n4d. 矩阵计算 (nfrac):")
vertexMatrix = cell.vertexMatrix
primaryVectors_INV = parameters.primaryVectors_INV
@benchmark $vertexMatrix * $primaryVectors_INV

println("\n5. 内存分配分析:")
println("5a. 单次Cascade的内存分配:")

# 创建新离子
ion3 = Atom(3, ionPosition, parameters)
SetVelocityDirection!(ion3, [0.,0.,-1.])
SetEnergy!(ion3, energy)
push!(simulator, ion3)

# 分析内存分配
allocs_before = Base.gc_alloc_count()
@time result = Cascade!(ion3, simulator)
allocs_after = Base.gc_alloc_count()
println("总内存分配: $(allocs_after - allocs_before) bytes")

println("\n6. 类型稳定性检查:")
println("检查关键函数的类型推断...")

# 检查LoadCell函数
println("6a. LoadCell:")
@code_warntype LoadCell(cell, simulator)

println("\n=== 性能优化建议 ===")
println("1. 查看 profile_results.txt 找出最耗时的函数")
println("2. 关注内存分配较多的地方") 
println("3. 检查类型不稳定的警告 (黄色/红色输出)")
println("4. 考虑预分配数组而不是使用push!")
println("5. 检查是否有不必要的重复计算")

println("\n=== 分析完成 ===") 