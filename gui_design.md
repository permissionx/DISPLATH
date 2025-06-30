# DISPLATH BCA Simulator - Modern GUI Design Document

## 设计理念
基于分析的代码结构，设计一个专业的BCA模拟器界面，注重：
- 科学计算的严谨性和专业性
- 参数输入的准确性和便利性  
- 模拟结果的清晰展示
- 现代化的视觉设计

## 界面布局架构

### 主布局：左右分栏 (70% + 30%)
```
┌─────────────────────────────────┬─────────────────┐
│           参数配置区域            │    结果监控区    │
│        (Parameter Setup)       │  (Results Panel)│
│                                │                 │
│  ┌─Material──Simulation──Ion─┐  │   ┌─Stats───┐   │
│  │         Tab Area          │  │   │         │   │
│  └─────────────────────────────┘  │   └─────────┘   │
│                                │   ┌─Progress─┐   │
│  ┌─Advanced Parameters──────┐   │   │         │   │
│  │  (Collapsible)           │   │   └─────────┘   │
│  └─────────────────────────────┘  │   ┌─Logs────┐   │
│                                │   │         │   │
│  ┌─Action Buttons───────────┐   │   │         │   │
│  │ [Validate] [Run] [Reset] │   │   │         │   │
│  └─────────────────────────────┘  │   └─────────┘   │
└─────────────────────────────────┴─────────────────┘
```

## 详细参数组织

### 1. Material Tab (材料配置)
```
Material Configuration
├── Preset Materials
│   └── Dropdown: [Custom, Silicon, MoS2, Graphene, ...]
├── Crystal Structure
│   ├── Primary Vectors (Å)
│   │   ├── a1: [4.260, 0.000, 0.000]
│   │   ├── a2: [0.000, 4.263, 0.000] 
│   │   └── a3: [0.000, 0.000, 6.700]
│   ├── Lattice Ranges (unit cells)
│   │   ├── X: [0] to [400]
│   │   ├── Y: [0] to [400]
│   │   └── Z: [2] to [1000]
│   └── Box Sizes (cells)
│       ├── nx: [400], ny: [400], nz: [1005]
│       └── Total atoms: ~400M (dynamic estimate)
├── Basis Atoms
│   ├── Atom 1: Position[0.000, 0.000, 0.000] Element[Si] ×
│   ├── Atom 2: Position[0.500, 0.500, 0.000] Element[Si] ×
│   ├── ... (expandable list)
│   └── [+ Add Atom]
└── Element Library
    ├── Si: Radius[2.20Å] Mass[28.09u] Z[14] DTE[35eV] ⚙️
    ├── Mo: Radius[1.40Å] Mass[95.96u] Z[42] DTE[60eV] ⚙️
    └── [+ Add Element]
```

### 2. Simulation Tab (模拟参数)
```
Simulation Parameters
├── Physical Parameters
│   ├── Temperature: [300.0] K
│   ├── Debye Temperature: [645.0] K
│   ├── Stop Energy: [20.0] eV
│   └── Vacancy Recover Distance: [0.0] Å
├── Computational Setup
│   ├── pMax: [1.45] Å
│   ├── Grid Vectors (Å)
│   │   ├── [6.0, 0.0, 0.0]
│   │   ├── [0.0, 6.0, 0.0]
│   │   └── [0.0, 0.0, 6.0]
│   ├── Dynamic Loading: ☑️ Enabled
│   ├── Cascades per Load: [10]
│   └── Random Seed: [43] [🎲 Random]
├── Output Settings
│   ├── Dump in Cascade: ☐ Disabled
│   ├── Output Directory: [./output/] [📁]
│   └── Log Level: [Info ▼]
└── Boundary Conditions
    ├── X Periodic: ☑️
    ├── Y Periodic: ☑️ 
    └── Z Periodic: ☐
```

### 3. Ion Tab (离子配置)
```
Ion Configuration
├── Ion Properties
│   ├── Ion Type: [Ar ▼] (from Element Library)
│   ├── Energy: [100000.0] eV
│   └── Energy Distribution: [Fixed ▼] [Gaussian] [Range]
├── Initial Conditions
│   ├── Position (Å)
│   │   ├── X: [Random in circle ▼] Radius[30.0]
│   │   ├── Y: [Random in circle ▼] Radius[30.0]
│   │   └── Z: [Fixed ▼] Value[2047.0]
│   └── Direction (normalized)
│       ├── θ (deg): [7.0] ├── φ (deg): [0.0]
│       └── Vector: [0.12, 0.00, -1.00]
├── Irradiation Protocol
│   ├── Number of Runs: [1000]
│   ├── Energy Steps: [10] from [10keV] to [100keV]
│   └── Flux: [2.0×10¹⁴] ions/cm²·s
└── Analysis Setup
    ├── Track Sputtering: ☑️
    ├── Track Defects: ☑️
    └── Track Range: ☑️
```

### 4. Advanced (可折叠高级参数)
```
Advanced Parameters (Collapsible)
├── θτ Data Configuration
│   ├── Repository Path: [/path/to/thetatau_repository/] 
│   ├── Energy Range: [-1.0:0.045:8.0]
│   ├── Impact Parameter Range: [0.0:0.01:5.0]
│   └── [🔄 Regenerate θτ Data]
├── Displacement Threshold Energy
│   ├── Mode: [Direct ▼] [Environment] [SOAP] [Custom]
│   ├── DTE File: [custom.dte] [📁]
│   └── Environment Cut: [3.0] Å
└── KMC Integration
    ├── Enable KMC: ☐
    ├── ν₀ Dictionary: [Edit] [📝]
    ├── Perfect Environment Index: [0]
    └── Irradiation Frequency: [0.0] Hz
```

## 右侧结果监控区 (30%)

### 1. Real-time Statistics
```
┌─ Simulation Progress ────────────┐
│ Run: 245/1000 (24.5%)           │
│ ████████░░░░░░░░░░░░░░ 24.5%     │
│ Current Energy: 78.5 keV        │
│ ETA: 15min 23s                  │
└─────────────────────────────────┘

┌─ Current Statistics ────────────┐
│ Total Collisions: 15,847       │
│ Displacements: 342             │
│ Vacancies: 156                 │
│ Interstitials: 186            │
│ Sputtered Atoms: 23            │
│ Cascade Depth: 847.3 Å        │
└─────────────────────────────────┘
```

### 2. Live Visualization
```
┌─ 3D Visualization ──────────────┐
│                                │
│    ● ● ● ○ ●    ← Atoms        │
│  ● ● ○ ● ● ●    ○ Vacancies    │
│    ● ● ● ● ●    → Ion path     │
│                                │
│ [Rotate] [Zoom] [Reset View]   │
└─────────────────────────────────┘
```

### 3. Real-time Plots
```
┌─ Energy vs Depth ───────────────┐
│ E│                             │
│  │ \                           │
│  │  \                          │
│  │   \___                      │
│  └─────────────────> Depth     │
└─────────────────────────────────┘
```

### 4. System Logs
```
┌─ System Logs ───────────────────┐
│ [14:32:15] Simulation started   │
│ [14:32:16] Loading 2.4M atoms   │
│ [14:32:18] θτ functions loaded  │
│ [14:32:19] Run 1/1000 started   │
│ [14:32:20] Cascade completed    │
│ [14:32:20] 23 vacancies created │
│ ⋮                              │
└─ [Clear] [Export] [Filter] ─────┘
```

## 现代化视觉设计

### 配色方案 (科学仪器风格)
- 主背景: 深蓝灰 (#1a1d29)
- 卡片背景: 中性蓝 (#252836) 
- 输入框: 浅蓝灰 (#333647)
- 主色调: 科技蓝 (#3b82f6)
- 强调色: 量子绿 (#10b981)
- 警告色: 能量橙 (#f59e0b)
- 文字: 冷白色系 (#e2e8f0, #cbd5e1, #94a3b8)

### 组件设计特色
1. **参数输入框**: 带单位后缀，科学计数法支持
2. **矩阵输入**: 3D网格可视化，实时验证
3. **元素库**: 周期表样式，元素属性面板
4. **进度条**: 多层显示（总体/当前/子任务）
5. **图表**: 实时更新的科学图表
6. **3D视窗**: 原子级可视化

### 交互特性
- 参数验证：实时检查物理合理性
- 智能建议：根据材料推荐参数
- 一键预设：常用材料配置模板
- 批量操作：能量扫描，参数优化
- 断点续传：长时间模拟的状态保存
- 结果对比：多次模拟结果叠加显示

## 响应式布局
- 大屏 (>1400px): 左右70/30分栏
- 中屏 (1024-1400px): 左右65/35分栏  
- 小屏 (<1024px): 上下堆叠布局
- 移动端: 标签页切换模式

这个设计体现了BCA模拟的专业性，同时保持了现代化的用户体验。 