# 🚀 Quick Start Guide for DISPLAΘ GUI

## Step-by-Step Testing Instructions

### 1. **Install Dependencies**
```bash
cd /beegfs/home/xuke/Researches/Irradiation_Li-Tianzhao/4.DISPLATH/DISPLATH/
julia gui/install_gui.jl
```

### 2. **Run Test Script** 
```bash
cd gui/
julia test_gui.jl
```
This will check all dependencies and core functionality.

### 3. **Launch GUI**
```bash
julia simple_gui.jl
```

## 🧪 **Testing Checklist**

### ✅ **Basic GUI Test**
1. **Window Opens**: GUI window should appear with tabs
2. **Material Selection**: Try changing between SiC, hBN, Graphene
3. **Parameter Input**: Enter values in different fields
4. **Tab Navigation**: Switch between Material, Simulation, Ion tabs

### ✅ **Quick Simulation Test**
**Use these minimal parameters for fast testing:**

**Material Tab:**
- Material: SiC  
- Lattice Constant: 4.36
- Box Size X: 5
- Box Size Y: 5  
- Box Size Z: 10

**Simulation Tab:**
- pMax: 4.0
- Vacancy Recover Distance: 4.0
- Temperature: 300.0
- Debye Temperature: 490.0
- Stop Energy: 10.0

**Ion Tab:**
- Ion Type: N
- Ion Energy: 1000.0 (lower energy for testing)
- Number of Cascades: 5 (very few for quick test)

### ✅ **Expected Results**
- Progress bar should advance
- Output log should show messages
- Status should update
- Results CSV file should be created
- Average vacancy count should be displayed

## 🔧 **Troubleshooting**

### **Problem: "Package Gtk not found"**
```bash
julia -e 'using Pkg; Pkg.add("Gtk")'
```

### **Problem: "JSON3 not found"**  
```bash
julia -e 'using Pkg; Pkg.add("JSON3")'
```

### **Problem: GUI doesn't open**
Check if you have a display:
```bash
echo $DISPLAY
```
If empty, you need X11 forwarding or a desktop environment.

### **Problem: "θτRepository not found"**
Check if the directory exists:
```bash
ls ../thetatau_repository/
```

### **Problem: Simulation crashes**
1. Try smaller box sizes (5x5x10)
2. Use fewer cascades (1-5)
3. Check output log for error messages

## 🎯 **Testing Scenarios**

### **Scenario 1: Basic Function Test**
- Open GUI ✓
- Change a parameter ✓  
- Save parameters ✓
- Load parameters ✓

### **Scenario 2: Quick Simulation**
- Set minimal parameters (above) ✓
- Run 1-2 cascades ✓
- Check results appear ✓

### **Scenario 3: Parameter Persistence**
- Set custom parameters ✓
- Save to JSON ✓
- Restart GUI ✓
- Load parameters ✓

## 📊 **What Success Looks Like**

### **GUI Launch Success:**
```
Starting DISPLAΘ GUI...
[GUI window opens with 3 tabs]
Status: Ready
```

### **Simulation Success:**
```
=== Starting DISPLAΘ Simulation ===
Parameters loaded: SiC with N ions
Creating simulator...
Running cascade 1/5
Running cascade 2/5
...
=== Simulation Complete ===
Average vacancies per cascade: 2.4
Results saved to: gui_results_2024-06-22T10:30:45.csv
```

## 🐛 **Common Issues & Solutions**

| Issue | Cause | Solution |
|-------|--------|----------|
| Gtk not found | Missing dependency | `Pkg.add("Gtk")` |
| GUI freezes | Too many cascades | Use 1-10 for testing |
| No display | Headless system | Use X11 forwarding |
| Simulation fails | Wrong parameters | Check θτRepository path |
| Memory error | Box too large | Start with 5x5x10 |

## 📞 **Getting Help**

If tests fail:
1. Check the test output for specific errors
2. Verify all dependencies are installed  
3. Try the minimal parameters first
4. Check file paths and permissions
5. Look at the output log for detailed error messages

## 🎉 **Next Steps After Successful Test**

1. **Try different materials** (when implemented)
2. **Experiment with parameters** 
3. **Run longer simulations**
4. **Analyze results in CSV files**
5. **Save/load parameter sets**

---

**Ready to test? Start with:**
```bash
cd gui/
julia test_gui.jl
```