# 开始计时
println("开始分子动力学模拟...")
overall_start = time()

home = ENV["ARCS_HOME"]
const BAlpha = 1.5
const BBeta = 0.44
const IS_DYNAMIC_LOAD = true

# 测量初始化时间
println("正在初始化...")
init_start = time()
include(home * "/src/DISPLATH.jl")
init_time = time() - init_start
println("初始化时间: $(round(init_time, digits=2)) 秒")

NI = 100000
flux = 5E14 * 1E-16
# Box and atoms 
a = 5.431
primaryVectors = [a 0.0 0.0; 0.0 a 0.0; 0.0 0.0 a]
boxL = ceil(Int, sqrt(NI/flux/primaryVectors[1,1]/primaryVectors[2,2]))
log_info("boxL: $boxL")
boxSizes = [boxL, boxL, 15005]
inputGridVectors = [5.0 0.0 0.0; 0.0 5.0 0.0; 0.0 0.0 5.0]
latticeRanges = [0 boxL; 0 boxL; 2 15000]
# Si diamond structure - conventional cell with 8 atoms
basis = [0.0 0.0 0.0;           # atom 1
         0.5 0.5 0.0;           # atom 2  
         0.5 0.0 0.5;           # atom 3
         0.0 0.5 0.5;           # atom 4
         0.25 0.25 0.25;        # atom 5
         0.75 0.75 0.25;        # atom 6
         0.75 0.25 0.75;        # atom 7
         0.25 0.75 0.75]        # atom 8
basisTypes = [1, 1, 1, 1, 1, 1, 1, 1]  # All 8 positions are Si atoms

# Parameters
pMax = a/3
vacancyRecoverDistance = 0.0
typeDict = Dict(
    1 => Element("Si", 20.0, 10.0),  # Si element
    2 => Element("B", 0.1, 0.1)     # B for ion bombardment    
)
seed = 43
const THREAD_RNG = [StableRNG(seed + t) for t in 1:Threads.nthreads()]

temperature = 300.0
DebyeTemperature = 519.0
isDumpInCascade = false
stopEnergy = 10.0
nCascadeEveryLoad = 1500
parameters = Parameters(primaryVectors, latticeRanges, basisTypes, basis,
                        pMax, vacancyRecoverDistance, typeDict; 
                        temperature=temperature, DebyeTemperature=DebyeTemperature, 
                        isDumpInCascade=isDumpInCascade, stopEnergy=stopEnergy, 
                        nCascadeEveryLoad=nCascadeEveryLoad, isAmorphous=false, amorphousLength=5.0, maxRSS=42)

# Run
simulator = Simulator(boxSizes, inputGridVectors, parameters)

phi = 30/180*π
thetas = [0, 2 ,4 ,7, 10]/180*π
n = 1
theta = thetas[n]
energy = 10_000.0
Restore!(simulator)

# 测量主循环时间
println("开始主循环 ($NI 次迭代)...")
loop_start = time()

cascade_times = Float64[]

for i in 1:NI
    # 进度报告
    if i % 4000 == 0 || i <= 10
        elapsed = time() - loop_start
        rate = i / elapsed
        remaining = (NI - i) / rate
        progress_percent = round(i / NI * 100, digits=1)
        println("进度: $i/$NI ($progress_percent%), 速率: $(round(rate, digits=2)) iter/s, 预计剩余: $(round(remaining, digits=2)) 秒")
    end
    
    offset = [0.0, 0.0, latticeRanges[3,2]*a-2]
    rp = RandomInSquare(boxSizes[1] * a, boxSizes[2] * a)
    ionPosition = Vector{Float64}(rp) + offset
    ion = Atom(2, ionPosition, parameters)
    v = [cos(phi)*sin(theta), sin(phi)*sin(theta), -cos(theta)]
    v = RandomlyDeviatedVector(v, 0.5/180.0*π)
    SetVelocityDirection!(ion, v)
    SetEnergy!(ion,energy)
    push!(simulator, ion)
    
    # 单独测量 Cascade! 时间
    cascade_start = time()
    Cascade!(ion, simulator)
    cascade_time = time() - cascade_start
    push!(cascade_times, cascade_time)
    
    z = ion.coordinate[3]
    @record  "R_p.csv" "$(a*latticeRanges[3,2]-z)" 
end

loop_time = time() - loop_start

# 输出计时结果
total_time = time() - overall_start

println("\n" * "="^50)
println("模拟统计报告")
println("="^50)
println("总运行时间: $(round(total_time, digits=2)) 秒")
println("初始化时间: $(round(init_time, digits=2)) 秒 ($(round(init_time/total_time*100, digits=1))%)")
println("主循环时间: $(round(loop_time, digits=2)) 秒 ($(round(loop_time/total_time*100, digits=1))%)")

if !isempty(cascade_times)
    using Statistics  # 确保统计函数可用
    
    println("\nCascade! 性能统计:")
    println("  平均每次迭代: $(round(loop_time/NI*1000, digits=2)) 毫秒")
    println("  最快 Cascade!: $(round(minimum(cascade_times)*1000, digits=2)) 毫秒")
    println("  最慢 Cascade!: $(round(maximum(cascade_times)*1000, digits=2)) 毫秒")
    println("  平均 Cascade!: $(round(mean(cascade_times)*1000, digits=2)) 毫秒")
    println("  Cascade! 时间标准差: $(round(std(cascade_times)*1000, digits=2)) 毫秒")
    
    # 计算性能分布
    sorted_times = sort(cascade_times)
    p50 = sorted_times[ceil(Int, length(sorted_times)*0.5)] * 1000
    p90 = sorted_times[ceil(Int, length(sorted_times)*0.9)] * 1000
    p99 = sorted_times[ceil(Int, length(sorted_times)*0.99)] * 1000
    println("  Cascade! 时间百分位数:")
    println("    中位数 (P50): $(round(p50, digits=2)) 毫秒")
    println("    P90: $(round(p90, digits=2)) 毫秒")
    println("    P99: $(round(p99, digits=2)) 毫秒")
end

println("\n最终输出...")
@dump  "defects.dump" [simulator.atoms; simulator.vacancies] ["vx", "vy", "vz", "e"]

final_time = time() - overall_start
println("完整模拟完成! 总耗时: $(round(final_time, digits=2)) 秒")