using LinearAlgebra
import Base: push!
using StableRNGs, Random, Base.Threads
using QuadGK
import Base: push!
import Base: delete!
using LinearAlgebra
using Interpolations   
using Dates
using ProgressMeter
using SparseArrays


#using PyCall
# @pyimport dscribe.descriptors as descriptors
#include("debug.jl")
#using .Recorder
include("types.jl")
include("elements.jl")
include("bca.jl")  # In namespace BCA: BCA->(QLoss, Constants)
using .BCA.ConstantFunctions
include("io.jl")
using .Output
include("geometry.jl")
include("dte.jl")
include("kmc.jl")
include("dynamics.jl")
include("dynamic_load.jl")
include("logging.jl")
include("utils.jl")


# modules


#simulator, atom_p, atom_t = test_collision()

# todo:
# Permutation at first in static load 

# 近似：
# 1. 多个碰撞一个
# 2. 间隙子位移能
# 3. 替代位位移能

# 未解决：
# 1. 入射粒子的位移能
# 2. 间隙子的坐标

# 增加内存监控与限制
# GUI
# 静态代码测试，因为改了cellIndex的类型为元组