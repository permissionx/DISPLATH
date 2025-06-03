# 快速性能检查脚本
using Profile

home = "/beegfs/home/xuke/Researches/Irradiation_Li-Tianzhao/4.DISPLATH/DISPLATH/"
include(home * "src/main.jl")

println("=== 快速性能检查 ===")

# 基本设置
a = 4.36
primaryVectors = [a 0.0 0.0; 0.0 a 0.0; 0.0 0.0 a]
boxSizes = [50, 50, 1400]
inputGridVectors = [a*2.1 0.0 0.0; 0.0 a*2.1 0.0; 0.0 0.0 a*2.1]
latticeRanges = [0 50; 0 50; 2 1200]
basis = [0.0 0.0 0.0; 0.25 0.25 0.25]
basisTypes = [1, 2]

typeDict = Dict(1 => Element("Si", 43.0, 22.0), 2 => Element("C", 40.0, 20.0), 3 => Element("N", 1.0, 1.0))
parameters = Parameters(primaryVectors, latticeRanges, basisTypes, basis, 
                       home * "thetatau_repository/", 4.0, 4.0, typeDict; 
                       temperature=300.0, DebyeTemperature=490.0)

println("1. 初始化...")
@time simulator = Simulator_dynamicLoad(boxSizes, inputGridVectors, parameters)
Save!(simulator)

println("\n2. 测试单个LoadCell函数:")
cell = simulator.cellGrid.cells[25, 25, 600]  # 中间位置的cell
cell.isLoaded = false
@time LoadCell(cell, simulator)
println("   生成原子数: $(length(cell.latticeAtoms))")

println("\n3. 测试5次LoadCell (不同位置):")
total_time = 0.0
total_atoms = 0
for i in 1:5
    test_cell = simulator.cellGrid.cells[20+i, 20+i, 500+i*10]
    test_cell.isLoaded = false
    time_result = @elapsed LoadCell(test_cell, simulator)
    atoms_count = length(test_cell.latticeAtoms)
    total_time += time_result
    total_atoms += atoms_count
    println("   Cell $(i): $(time_result*1000)ms, $(atoms_count) atoms")
end
println("   平均: $(total_time/5*1000)ms/cell, $(total_atoms/5) atoms/cell")

println("\n4. Profile单次Cascade:")
# 创建离子
ion = Atom(3, [110.0, 110.0, 5240.0], parameters)
SetVelocityDirection!(ion, [0.,0.,-1.])
SetEnergy!(ion, 10000.0)
push!(simulator, ion)

Profile.clear()
@profile Cascade_dynamicLoad!(ion, simulator)

println("\n=== Top函数 (>1% 时间) ===")
Profile.print(maxdepth=10, mincount=100)

# 保存详细结果
open("quick_profile.txt", "w") do io
    println(io, "=== 快速Profile结果 ===")
    Profile.print(io, maxdepth=15, mincount=50)
end

println("\n快速检查完成！详细结果保存在 quick_profile.txt") 