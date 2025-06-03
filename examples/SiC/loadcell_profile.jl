# LoadCell函数专门性能分析
using Profile
using BenchmarkTools

home = "/beegfs/home/xuke/Researches/Irradiation_Li-Tianzhao/4.DISPLATH/DISPLATH/"
include(home * "src/main.jl")

println("=== LoadCell函数性能分析 ===")

# 初始化
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

simulator = Simulator_dynamicLoad(boxSizes, inputGridVectors, parameters)

# 选择测试cell
cell = simulator.cellGrid.cells[25, 25, 600]

println("1. LoadCell整体性能:")
cell.isLoaded = false
@benchmark LoadCell($cell, $simulator) setup=(cell.isLoaded = false)

println("\n2. LoadCell各步骤详细分析:")

# 手动执行LoadCell中的每个步骤来测量
cell.isLoaded = false
println("   参数提取:")
@time begin
    parameters_local = simulator.parameters
    primaryVectors_INV = parameters_local.primaryVectors_INV
    primaryVectors_local = parameters_local.primaryVectors
    latticeRanges_local = parameters_local.latticeRanges
    basisTypes_local = parameters_local.basisTypes
    basis_local = parameters_local.basis
    vertexMatrix = cell.vertexMatrix
end

println("   矩阵乘法 (nfrac):")
@time nfrac = vertexMatrix * primaryVectors_INV

println("   计算nmin/nmax:")
@time begin
    nmin = [floor(Int, minimum(nfrac[:, 1]) - 1), 
            floor(Int, minimum(nfrac[:, 2]) - 1), 
            floor(Int, minimum(nfrac[:, 3]) - 1)]
    nmax = [ceil(Int, maximum(nfrac[:, 1]) + 1),
            ceil(Int, maximum(nfrac[:, 2]) + 1),
            ceil(Int, maximum(nfrac[:, 3]) + 1)]
end

println("   计算范围:")
@time begin
    n1 = max(nmin[1],latticeRanges_local[1,1]):min(nmax[1],latticeRanges_local[1,2])
    n2 = max(nmin[2],latticeRanges_local[2,1]):min(nmax[2],latticeRanges_local[2,2])
    n3 = max(nmin[3],latticeRanges_local[3,1]):min(nmax[3],latticeRanges_local[3,2])
    ranges = cell.ranges
    x_min, x_max = ranges[1, 1], ranges[1, 2]
    y_min, y_max = ranges[2, 1], ranges[2, 2] 
    z_min, z_max = ranges[3, 1], ranges[3, 2]
end

println("   原子生成循环:")
@time begin
    atoms = Vector{Atom}()
    idx = -1
    for x in n1 
        for y in n2
            for z in n3
                for i in basisTypes_local
                    coordinate = primaryVectors_local' * (Float64[x, y, z] + basis_local[i, :])
                    
                    if (coordinate[1] >= x_min && coordinate[1] <= x_max &&
                        coordinate[2] >= y_min && coordinate[2] <= y_max &&
                        coordinate[3] >= z_min && coordinate[3] <= z_max)
                        
                        atom = Atom(basisTypes_local[i], coordinate, parameters_local)
                        atom.cellIndex = cell.index
                        atom.index = idx
                        idx -= 1
                        push!(atoms, atom)
                    end
                end
            end
        end
    end
end

println("   生成原子数: $(length(atoms))")

println("\n3. 各部分微基准测试:")

println("   3a. 矩阵乘法:")
@benchmark $vertexMatrix * $primaryVectors_INV

println("\n   3b. minimum/maximum计算:")
@benchmark begin
    [floor(Int, minimum($nfrac[:, 1]) - 1), 
     floor(Int, minimum($nfrac[:, 2]) - 1), 
     floor(Int, minimum($nfrac[:, 3]) - 1)]
end

println("\n   3c. 单个原子创建:")
test_coord = [100.0, 100.0, 5000.0]
@benchmark Atom(1, $test_coord, $parameters)

println("\n   3d. 坐标变换:")
test_fractional = [20.0, 20.0, 500.0]
@benchmark $primaryVectors_local' * $test_fractional

println("\n4. 内存分配分析:")
function analyze_loadcell_memory(cell, simulator)
    cell.isLoaded = false
    LoadCell(cell, simulator)
    return length(cell.latticeAtoms)
end

cell.isLoaded = false
@time atom_count = analyze_loadcell_memory(cell, simulator)
println("   生成原子数: $atom_count")

println("\n5. 优化建议:")
println("   - 如果矩阵乘法慢：考虑预计算或缓存")
println("   - 如果原子创建慢：考虑预分配Vector并批量创建")
println("   - 如果循环慢：考虑向量化操作")
println("   - 如果边界检查慢：考虑优化判断逻辑")

println("\n=== 分析完成 ===") 