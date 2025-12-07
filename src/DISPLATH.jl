"""
DISPLATH - 基于Julia的二元碰撞近似(BCA)模拟器

用于模拟材料中的离子辐照过程，包括碰撞级联、缺陷统计和动力学蒙特卡洛模拟。

主要功能：
- 二元碰撞近似(BCA)物理模型
- 碰撞级联模拟
- 缺陷统计与分析
- 动态加载模式（节省内存）
- 动力学蒙特卡洛(KMC)模拟

使用示例：
```julia
using DISPLATH
# 创建参数和模拟器
parameters = Parameters(...)
simulator = Simulator(...)
# 运行模拟
Cascade!(ion, simulator)
```
"""
module DISPLATH

# 外部依赖
using LinearAlgebra
using StaticArrays
import Base: push!, delete!
using StableRNGs, Random, Base.Threads
using QuadGK
using Interpolations   
using Dates
using ProgressMeter
using Distributions

# 包含核心文件（按依赖顺序）
include("debug.jl")
include("logging.jl")  # Logging 模块 - 提前包含以便在 types.jl 中使用异常类型
using .Logging
include("types.jl")
include("elements.jl")
include("bca.jl")  # BCA 模块
using .BCA.ConstantFunctions
include("io.jl")    # Output 模块
using .Output
include("geometry.jl")
include("dte.jl")
include("kmc.jl")
include("dynamics.jl")
include("dynamic_load.jl")
include("utils.jl")
include("interface.jl")

# ========== 公共 API 导出 ==========
# 核心类型
export Atom, Simulator, Parameters, Cell, Grid, Box, Material
export LatticePoint, Element, ConstantsByType
export MigrationEvent  # 如果存在原子扩散模块

# 核心函数
export Cascade!, Restore!, Save!
export DefectStatics, CountVacancies
export GetCell, WhichCell, ComputeDistance
export CreateGrid, CreateBoxByPrimaryVectors
export Irradiation  # 通用辐照函数

# 统一接口函数（自动选择静态/动态加载实现）
export ShotTarget, Collision!, DumpInCascade, Dump
export LeaveLatticePoint!, Stop!

# 原子操作
export SetEnergy!, SetVelocityDirection!
export DisplaceAtom!, SetOnLatticePoint!

# 材料与模拟器创建
export Material, Simulator

# 工具函数
export RandomPointInCircle, RandomInSquare, RandomVectorInUnitSphere
export RandomlyDeviatedVector, rotation_matrix_from_vectors

# 日志函数（从 logging.jl）
export log_info, log_success, log_warning, log_error, log_debug, log_section, log_separator

# 错误处理（从 logging.jl）
export DISPLATHError, InvalidParameterError, SimulationError, GeometryError

# I/O 宏（从 Output 模块）
export @dump, @record

# BCA 模块导出
export CollisionParams, Q_nl

end # module DISPLATH


# modules


#simulator, atom_p, atom_t = test_collision()

# todo:
# Permutation at first in static load 

# TODO: 
# 1. Multiple collisions with one target
# 2. Interstitial displacement energy  
# 3. Substitutional displacement energy

# Unresolved:
# 1. Incident particle displacement energy
# 2. Interstitial coordinates

# todo:
# Features:
# Added memory monitoring and limits: need to restore all atoms when clearing cells ... ok
# Dynamic load system can be set to infinite size
# Static code testing, changed cellIndex type to tuple ... ok  
# Multi-component simulation through simulator composition

