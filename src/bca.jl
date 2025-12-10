#=数值方法分析
1. 积分技术
自适应积分：处理近奇异被积函数
截断策略：无穷积分的有穷近似
精度控制：rtol=1e-8的相对误差容限
2. 插值优化
对数网格：在能量和碰撞参数空间使用对数采样
预计算：避免实时数值积分
精度平衡：计算速度与物理准确性的权衡
3. 数值稳定性
边界处理：小p值的特殊处理
溢出防护：能量损失的合理限制
收敛保证：二分法的严格收敛条件
物理模型总结
bca.jl实现了完整的经典二体碰撞物理：
势函数：ZBL Universal势
运动学：相对论前的经典力学
能量损失：核阻止 + 电子阻止
数值方法：自适应积分 + 插值优化
在DISPLATH架构中的地位
这个文件是DISPLATH的物理引擎核心：
碰撞检测：ComputeP!函数
运动学计算：CollisionParams函数
能量沉积：E_t和QLoss模块
缺陷产生：通过能量转移判断
这种设计使得DISPLATH能够在保持物理准确性的同时，实现高效的辐射损伤模拟。=#
module BCA
using LinearAlgebra
using QuadGK
using Base.MathConstants
using Main: ConstantsByType

export CollisionParams, Q_nl

const qe_squared = Float64(14.399764) # square of element charge, unit: eV*angstrom

#通用势函数 - Universal Potential
#= 物理意义：Ziegler-Biersack-Littmark (ZBL) 通用势的屏蔽函数

理论基础：基于量子力学计算的原子间势经验公式

数学形式：4项指数衰减函数的线性组合

适用范围：适用于所有原子组合，能量范围1eV-1MeV =#
function Φ(x::Float64)
    A = Vector{Float64}([0.1818, 0.5099, 0.2802, 0.02817])
    B = Vector{Float64}([3.2, 0.9423, 0.4028, 0.2016])
    result = 0.0
    for i in 1:4
        result += A[i] * exp(-B[i] * x)
    end
    return result
end

#原子间势函数
#=Z_p, Z_t：入射粒子和靶原子的原子序数

a_U：Universal势的屏蔽长度

Φ：无量纲屏蔽函数=#
function V(r::Float64, type_p::Int64, type_t::Int64, constantsByType::ConstantsByType)
    return constantsByType.V_upterm[[type_p, type_t]] / r * Φ(r / constantsByType.a_U[[type_p, type_t]])  
end

#约化能量计算
#=物理意义：质心系中的约化能量

公式：E_r = [m_t/(m_p + m_t)] × E_p

用途：将实验室系能量转换为质心系能量=#
function E_r(energy_p::Float64, mass_p::Float64, mass_t::Float64)
    return mass_t * energy_p / (mass_p + mass_t)
end

#径向运动方程
#=物理意义：经典二体问题的径向速度分量

方程：g(r) = √[1 - p²/r² - V(r)/E_r]

来源：从能量和角动量守恒推导得出=#
function g(r::Float64, p_squared::Float64, E_r::Float64, type_p::Int64, type_t::Int64, constantsByType::ConstantsByType)
    value = 1 - p_squared / (r * r) - V(r, type_p, type_t, constantsByType) / E_r
    return value > 0.0 ? sqrt(value) : 0.0
end

#转折点搜索算法
#=算法分析：

二分查找：在[guessStart, 100]范围内搜索转折点

收敛条件：1e-14 Å的精度要求

物理意义：找到g(r) = 0的点，即经典轨道的最接近距离=#
function FindTurningPoint(p_squared::Float64, E_r::Float64,type_p::Int64, type_t::Int64,  guessStart::Float64, constantsByType::ConstantsByType)
    rLeft = guessStart
    rRight = 100
    while (rRight - rLeft) > 1e-14
        rMid = (rLeft + rRight) / 2
        g_squared = 1 - p_squared / (rMid * rMid) - V(rMid, type_p, type_t, constantsByType) / E_r
        if g_squared < 0
            rLeft = rMid
        elseif g_squared > 0
            rRight = rMid
        else
            return rMid+1e-15
        end
    end
    return rRight
end

#散射角积分计算
#=物理理论：经典散射角计算公式
θ = π - 2p ∫[r_min→∞] dr / [r² × g(r)]
数值方法：

自适应高斯-克朗罗德积分：quadgk函数处理奇异积分

截断处理：积分到rStart+10000，然后加渐近项修正=#
function Integrate_g_θ(p_squared::Float64, E_r::Float64, type_p::Int64, type_t::Int64, rStart::Float64, constantsByType::ConstantsByType)
    result = quadgk(r -> 1 / (r * r * g(r, p_squared, E_r, type_p, type_t, constantsByType)),
        rStart, rStart+10000,
        rtol=1e-8)[1]
    return result+1/(rStart+10000)
end


function θ(p::Float64, p_squared::Float64, E_r::Float64, type_p::Int64, type_t::Int64, rStart::Float64, constantsByType::ConstantsByType)
    return π - 2 * p * Integrate_g_θ(p_squared, E_r, type_p, type_t, rStart, constantsByType)
end

##飞行时间计算
#=物理意义：碰撞过程中的时间延迟

数学处理：减去自由粒子项以处理积分发散

用途：计算碰撞点的精确位置=#
function Integrate_g_τ(p_squared::Float64, E_r::Float64, type_p::Int64, type_t::Int64, rStart::Float64, constantsByType::ConstantsByType)
    result = quadgk(r -> 1 / g(r, p_squared, E_r, type_p, type_t, constantsByType) - 1 / sqrt(1 - p_squared / (r * r)),
        rStart, rStart+1000,
        rtol=1e-8)[1]
    return result
end


function τ(p_squared::Float64, type_p::Int64, type_t::Int64, E_r::Float64, rStart::Float64, constantsByType::ConstantsByType)
    return sqrt(rStart * rStart - p_squared) - Integrate_g_τ(p_squared, E_r, type_p, type_t, rStart, constantsByType)
end

#碰撞后运动学计算
#=物理理论：实验室系中的散射角

tanφ：入射粒子的散射角正切

tanψ：反冲靶原子的反冲角正切

f：能量损失因子=#
function tanφ(mass_p::Float64, mass_t::Float64, θ::Float64, f::Float64)
    return mass_t * sin(θ) * f / (mass_p + mass_t * cos(θ) * f)
end

function tanψ(θ::Float64, f::Float64)
    return sin(θ) * f / (1 - cos(θ) * f)
end
#能量转移计算 物理公式：二体弹性碰撞的能量转移 最大能量转移：当θ=π时，E_t,max = [4 m_p m_t / (m_p + m_t)²] × E_p
function E_t(mass_p::Float64, mass_t::Float64, energy_p::Float64, θ::Float64)
    return 4 * mass_p * mass_t / (mass_p + mass_t)^2 * energy_p * sin(θ / 2)^2
end

function E_p(energy_p::Float64, E_t::Float64, Q::Float64)
    return energy_p - E_t - Q
end

#位移计算 
#=物理意义：碰撞点相对于原子位置的位移

x_p：入射粒子在速度方向的位移

x_t：靶原子在速度方向的位移=#
function x_p(mass_p::Float64,mass_t::Float64, p::Float64, θ::Float64, τ::Float64)
    return (2 * mass_p * τ + (mass_t - mass_p) * p * tan(θ / 2)) / (mass_p + mass_t)
end

function x_t(p::Float64, θ::Float64, x_p::Float64)
    return p * tan(θ / 2) - x_p
end

function θτ(energy_p::Float64, mass_p::Float64, mass_t::Float64, type_p::Int64, type_t::Int64,
    p::Float64, constantsByType::ConstantsByType)
    p_squared = p * p
    E_r_v = E_r(energy_p, mass_p, mass_t)
    rStart = FindTurningPoint(p_squared, E_r_v, type_p, type_t, p, constantsByType)
    θ_v = θ(p, p_squared, E_r_v, type_p, type_t, rStart, constantsByType)
    τ_v = τ(p_squared, type_p, type_t, E_r_v, rStart, constantsByType)
    return θ_v, τ_v
end

#主碰撞参数函数
#=关键优化：
插值使用：对于p > 1e-10使用预计算的插值函数
直接计算：对于小p值直接数值积分
能量损失：包含局域电子阻止本领效应=#
function CollisionParams(energy_p::Float64, mass_p::Float64, mass_t::Float64, type_p::Int64, type_t::Int64,
                         p::Float64, constantsByType::ConstantsByType,
                         θFunction::Function, τFunction::Function)
    E_r_v = E_r(energy_p, mass_p, mass_t)
    Q_loc_v = QLoss.Q_loc(energy_p, type_p, type_t, E_r_v, p, constantsByType)
    Q_loc_v = min(Q_loc_v, (1 - 1E-6) * E_r_v)
    f = sqrt(1 - Q_loc_v / E_r_v)
    #f = 1.0
    #rStart = FindTurningPoint(p_squared, E_r_v, type_p, type_t, p, constantsByType)
    #θ_v = θ(p, p_squared, E_r_v, type_p, type_t, rStart, constantsByType)
    #τ_v = τ(p_squared, type_p, type_t, E_r_v, rStart, constantsByType)
    #@show θ_v, τ_v
    if p > 1E-10
        ePPower = log10(energy_p)
        pPower = log10(p)
        θ_v = θFunction(ePPower, pPower)
        τ_v = τFunction(ePPower, pPower)
    else
        θ_v, τ_v = θτ(energy_p, mass_p, mass_t, type_p, type_t, p, constantsByType)
    end
    #@show θ_v, τ_v
    tanφ_v = tanφ(mass_p, mass_t, θ_v, f)
    tanψ_v = tanψ(θ_v, f)
    E_t_v = E_t(mass_p, mass_t, energy_p, θ_v)
    x_p_v = x_p(mass_p, mass_t, p, θ_v, τ_v)
    x_t_v = x_t(p, θ_v, x_p_v)
    return tanφ_v, tanψ_v,  E_t_v, x_p_v, x_t_v, Q_loc_v
end

#QLoss模块 - 电子能量损失
#=物理模型：Lindhard-Scharff电子阻止公式
理论基础：基于自由电子气模型的量子力学计算
能量依赖：覆盖低能核阻止到高能电子阻止的过渡=#
function Q_nl(energy_p::Float64, mass_p::Float64, mass_t::Float64, type_p::Int64, type_t::Int64,
                         pL::Float64, N::Float64, constantsByType::ConstantsByType)
    E_r_v = E_r(energy_p, mass_p, mass_t)
    Q_nl_v = QLoss.Q_nl(energy_p, type_p, type_t, E_r_v, pL, N, constantsByType)
    return Q_nl_v
end


module QLoss
using Main: ConstantsByType
# 
# constants for different types:
function S_e(energy_p::Float64, type_p::Int64, type_t::Int64, constantsByType::ConstantsByType)
    δ_1_2 = Float64(0.7125) # δ^(1/2)
    δ_1_2_ = Float64(-0.7125) # -δ^(1/2)
    δ_1 = Float64(0.7017543859649122) # 1/δ
    A = energy_p / constantsByType.E_m[type_p]
    termLeftDown = (A / log(A + 1 / A + ℯ - 2))^δ_1_2 
    termRightDown = A ^ δ_1_2_
    termUp = constantsByType.S_e_upTerm[[type_p, type_t]]
    return termUp / (termLeftDown + termRightDown)^δ_1
end


function x_nl(type_p::Int64, type_t::Int64, E_r::Float64, constantsByType::ConstantsByType)
    # termOthers 
    termE_r = E_r^0.075
    termOthers = constantsByType.x_nl[[type_p, type_t]]
    return termE_r * termOthers
end 

function x_loc(x_nl::Float64)
    return 1 - x_nl
end

#非局域能量损失
#=物理意义：沿离子路径积分的非局域电子能量损失

pL：飞行路径长度

N：局部原子密度

S_e：电子阻止本领=#
function Q_nl_f(type_p::Int64, type_t::Int64, S_e::Float64, x_nl::Float64, x_loc::Float64, pL::Float64, N::Float64, constantsByType::ConstantsByType)
    # constant:  pMax (half of lattice constant)
    termRight = x_nl + x_loc * constantsByType.Q_nl[[type_p, type_t]]
    return S_e * N * termRight * pL 
end 

#局域能量损失
#=物理模型：与碰撞参数相关的局域电子激发

指数衰减：exp(-p/a)表示随距离的衰减

特征长度：a为电子云的特征尺寸=#
function Q_loc_f(type_p::Int64, type_t::Int64, S_e::Float64, x_loc::Float64, p::Float64, constantsByType::ConstantsByType)
    termUp = x_loc * S_e * exp(-p / constantsByType.a[[type_p, type_t]])
    termDown = constantsByType.Q_loc[[type_p, type_t]]
    return termUp * termDown
end


function Q_loc(energy_p::Float64, type_p::Int64, type_t::Int64, E_r::Float64, p::Float64, constantsByType::ConstantsByType)
    S_e_v = S_e(energy_p, type_p, type_t, constantsByType)
    x_nl_v = x_nl(type_p, type_t, E_r, constantsByType)
    x_loc_v = x_loc(x_nl_v)
    Q_loc_v = Q_loc_f(type_p, type_t, S_e_v, x_loc_v, p, constantsByType)
    # println("S_e: ", S_e_v, "\n", "x_nl: ", x_nl_v, "\n", "x_loc: ", x_loc_v, "\n", "Q_nl: ", Q_nl_v, "\n", "Q_loc: ", Q_loc_v)
    return Q_loc_v 
end 

function Q_nl(energy_p::Float64, type_p::Int64, type_t::Int64, E_r::Float64, pL::Float64, N::Float64, constantsByType::ConstantsByType)
    S_e_v = S_e(energy_p, type_p, type_t, constantsByType)
    x_nl_v = x_nl(type_p, type_t, E_r, constantsByType)
    x_loc_v = x_loc(x_nl_v)
    Q_nl_v = Q_nl_f(type_p, type_t, S_e_v, x_nl_v, x_loc_v, pL, N, constantsByType)
    return Q_nl_v
end 

end

#物理常数
module ConstantFunctions
using Base.MathConstants
using ..BCA: qe_squared

#势函数参数计算
function V_upterm(Z_p::Float64, Z_t::Float64)
    return Z_p * Z_t * qe_squared
end

function a_U(Z_p::Float64, Z_t::Float64) #ZBL屏蔽长度公式，其中a_0 = 0.529 Å为玻尔半径
    a_0 = 0.529177210903::Float64  # Bohr radius, unit: angstrom
    return 0.8854 * a_0 / (Z_p^0.23 + Z_t^0.23)
end

#=物理意义：Lindhard特征能量

用于电子阻止本领的能量标度

与原子电子云密度相关=#
function E_m(Z_p::Float64, mass_p::Float64)#特征能量计算
    v_b_squared = Float64(49708.0)  # square of Bohr velocity, unit: eV/u need to check carefully
    return 0.5 * mass_p * v_b_squared * Z_p^(4/3)
end

function S_e_upTerm(type_p::Int64,Z_p::Float64, Z_t::Float64, m_p::Float64, α_p::Float64)
    # In initializations
    k_LS = 1.212 * Z_p^(7/6) * Z_t / ((Z_p^(2/3) + Z_t^(2/3))^(3/2) * m_p^(1/2))
    return α_p * k_LS * sqrt(E_m(Z_p, m_p))
end

function x_nl(type_p::Int64, Z_p::Float64, Z_t::Float64, β_p::Float64)
    return β_p * (a_U(Z_p, Z_t) / (Z_p * Z_t * qe_squared))^0.075
end

function a(Z_p::Float64, Z_t::Float64)
    return 1.45 * a_U(Z_p, Z_t) / 0.3 / Z_p^0.4
end

function Q_nl(Z_p::Float64, Z_t::Float64, pMax::Float64)
    return (1 + pMax / a(Z_p, Z_t)) * exp(-pMax / a(Z_p, Z_t))
end

function Q_loc(Z_p::Float64, Z_t::Float64)
    return 1 / (2 * π * a(Z_p, Z_t)^2)
end

end


end
