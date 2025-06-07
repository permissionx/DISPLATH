home = "/beegfs/home/xuke/Researches/Irradiation_Li-Tianzhao/4.DISPLATH/DISPLATH/"
include(home * "src/main.jl")
using Random
using Profile
using Dates

# Box and atoms 
a = 4.36
primaryVectors = [a 0.0 0.0; 0.0 a 0.0; 0.0 0.0 a]
boxSizes = [800, 800, 2002]
inputGridVectors = [a*4.1 0.0 0.0; 0.0 a*4.1 0.0; 0.0 0.0 a*4.1]
latticeRanges = [0 800; 0 800; 2 2000]
basis = [0.0 0.0 0.0; 0.25 0.25 0.25]
basisTypes = [1, 2]


# Parameters
θτRepository = home * "thetatau_repository/"
pMax = 4.0
vacancyRecoverDistance = 10.0
typeDict = Dict(
    1 => Element("Si", 20.0, 10.0),  
    2 => Element("C", 35.0, 17.5),
    3 => Element("N", 1.0, 1.0)    
)
temperature = 300.0
DebyeTemperature = 300.0
isDumpInCascade = true
isDynamicLoad = true
stopEnergy = 20.0
parameters = Parameters(primaryVectors, latticeRanges, basisTypes, basis,
                        θτRepository, pMax, vacancyRecoverDistance, typeDict; 
                        temperature=temperature, DebyeTemperature=DebyeTemperature, 
                        isDumpInCascade=isDumpInCascade, isDynamicLoad=isDynamicLoad,
                        stopEnergy=stopEnergy)

# Run
simulator = Simulator_dynamicLoad(boxSizes, inputGridVectors, parameters)
#cell = simulator.cellGrid.cells[10, 10, 10]
#@time LoadCell(cell, simulator)

Random.seed!(43)
energy = 28000.0
#ionPosition = RandomPointInCircle(20.0) + [110.0, 110.0, 5240.0] 
#ion = Atom(3, ionPosition, parameters)
#SetVelocityDirection!(ion, [0.,0.,-1.])
#SetEnergy!(ion,energy)
#push!(simulator, ion)
#Cascade_dynamicLoad!(ion, simulator)

#Profile.clear()
println("Running 1000 ions...")

# 初始内存状态
initial_memory = Sys.maxrss()
println("Initial memory usage: $(round(initial_memory/1024/1024, digits=2)) MB")

# 创建内存监控日志文件
memory_log_file = "memory_usage.txt"
open(memory_log_file, "w") do io
    println(io, "=== Memory Usage Monitor ===")
    println(io, "Start time: $(Dates.now())")
    println(io, "Initial memory usage: $(round(initial_memory/1024/1024, digits=2)) MB")
    println(io, "Memory check interval: every 50 iterations")
    println(io, "=" ^ 50)
    println(io, "")
end

if simulator.parameters.isDumpInCascade
    Dump_dynamicLoad(simulator, "Cascade.dump", 0, false)
end

memory_check_interval = 50  # 每50次检查一次内存

@showprogress for i in 1:1000
    ionPosition = RandomPointInCircle(400.0) + [1740.0, 1740.0, 8726.0] 
    ion = Atom(3, ionPosition, parameters)
    SetVelocityDirection!(ion, [0.0,0.,-1.])
    SetEnergy!(ion,energy)
    push!(simulator, ion)
    Cascade_dynamicLoad!(ion, simulator)
    
    # 每隔指定次数检查内存使用量
    if i % memory_check_interval == 0
        current_memory = Sys.maxrss()
        gc_bytes = Base.gc_live_bytes()
        memory_increase = current_memory - initial_memory
        
        # 输出到控制台（简化版）
        println("Iteration $i: Memory $(round(current_memory/1024/1024, digits=2)) MB (+$(round(memory_increase/1024/1024, digits=2)) MB)")
        
        # 详细信息输出到文件
        open(memory_log_file, "a") do io
            println(io, "Iteration $i ($(Dates.now())):")
            println(io, "  Current memory: $(round(current_memory/1024/1024, digits=2)) MB")
            println(io, "  Memory increase: $(round(memory_increase/1024/1024, digits=2)) MB")
            println(io, "  GC live bytes: $(round(gc_bytes/1024/1024, digits=2)) MB")
            println(io, "  Memory efficiency: $(round(gc_bytes/current_memory*100, digits=1))%")
            println(io, "")
        end
        
        # 可选：强制垃圾回收
        # GC.gc()
    end
end

# 程序结束时记录最终状态
final_memory = Sys.maxrss()
open(memory_log_file, "a") do io
    println(io, "=" ^ 50)
    println(io, "Final memory usage: $(round(final_memory/1024/1024, digits=2)) MB")
    println(io, "Total memory increase: $(round((final_memory-initial_memory)/1024/1024, digits=2)) MB")
    println(io, "End time: $(Dates.now())")
end

println("Memory monitoring completed. Check '$memory_log_file' for detailed results.")
#Profile.print(maxdepth=12, mincount=100)

#println("\n保存详细Profile结果到文件...")
#open("cascade_profile.txt", "w") do io
#    println(io, "=== Cascade_dynamicLoad Profile结果 ===")
#    println(io, "")
#    Profile.print(io, maxdepth=20, mincount=50)
#end







