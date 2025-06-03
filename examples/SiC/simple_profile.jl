# 简单的Cascade_dynamicLoad性能分析
using Profile

home = "/beegfs/home/xuke/Researches/Irradiation_Li-Tianzhao/4.DISPLATH/DISPLATH/"
include(home * "src/main.jl")

println("=== Cascade_dynamicLoad 性能分析 ===")

# 基本设置 (复制自run.jl)
a = 4.36
primaryVectors = [a 0.0 0.0; 0.0 a 0.0; 0.0 0.0 a]
boxSizes = [50, 50, 1400]
inputGridVectors = [a*2.1 0.0 0.0; 0.0 a*2.1 0.0; 0.0 0.0 a*2.1]
latticeRanges = [0 50; 0 50; 2 1200]
basis = [0.0 0.0 0.0; 0.25 0.25 0.25]
basisTypes = [1, 2]

typeDict = Dict(
    1 => Element("Si", 43.0, 22.0),  
    2 => Element("C", 40.0, 20.0),
    3 => Element("N", 1.0, 1.0)    
)

parameters = Parameters(primaryVectors, latticeRanges, basisTypes, basis, 
                       home * "thetatau_repository/", 4.0, 4.0, typeDict; 
                       temperature=300.0, DebyeTemperature=490.0)

println("1. 初始化模拟器...")
@time simulator = Simulator_dynamicLoad(boxSizes, inputGridVectors, parameters)

println("\n2. 运行3个离子 (计时):")
energy = 10000.0
total_time = 0.0

for i in 1:3
    println("   离子 $i:")
    ionPosition = [110.0, 110.0, 5240.0] 
    ion = Atom(3, ionPosition, parameters)
    SetVelocityDirection!(ion, [0.,0.,-1.])
    SetEnergy!(ion, energy)
    push!(simulator, ion)
    
    global total_time
    time_result = @elapsed Cascade_dynamicLoad!(ion, simulator)
    total_time += time_result
    println("     耗时: $(time_result) 秒")
end

println("   平均耗时: $(total_time/3) 秒/离子")

println("\n3. Profile分析 (运行1个离子):")
# 创建新离子进行profile
ionPosition = [110.0, 110.0, 5240.0] 
ion_profile = Atom(3, ionPosition, parameters)
SetVelocityDirection!(ion_profile, [0.,0.,-1.])
SetEnergy!(ion_profile, energy)
push!(simulator, ion_profile)

# Profile分析
Profile.clear()
@profile Cascade_dynamicLoad!(ion_profile, simulator)

println("\n=== Profile结果 (主要瓶颈) ===")
Profile.print(maxdepth=12, mincount=100)

# 保存详细结果
println("\n保存详细Profile结果到文件...")
open("cascade_profile.txt", "w") do io
    println(io, "=== Cascade_dynamicLoad Profile结果 ===")
    println(io, "运行了$(simulator.numberOfAtoms)个原子")
    println(io, "平均耗时: $(total_time/3) 秒/离子")
    println(io, "")
    Profile.print(io, maxdepth=20, mincount=50)
end

println("\n=== 内存统计 ===")
println("GC时间: $(Base.gc_time_ns()/1e9) 秒")
println("总原子数: $(simulator.numberOfAtoms)")

println("\n=== 分析完成 ===")
println("详细结果已保存到 cascade_profile.txt")
println("查看主要瓶颈请关注上面Profile结果中耗时最多的函数") 