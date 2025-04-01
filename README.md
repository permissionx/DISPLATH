# DISPLAΘ

> **Note**: This README was generated by GPT. A comprehensive manual documentation will be written later.

DISPLAΘ is a Julia-based Binary Collision Approximation (BCA) simulator designed to model ion irradiation processes in materials. The name DISPLAΘ sounds like "displace" and stands for "DISPLAΘ Is the Second Program Laveraging Accurate Θ" (where I2DM [1] was the first one). Through the simulation of collision cascades, it helps predict damage such as vacancies produced by different ions at various energies. ⚡

The name DISPLAΘ is a clever acronym that reflects the program's core functionality:
- **DISPLA**: Represents the program's ability to model displacement and damage in materials
- **Θ**: The Greek letter theta, symbolizing the accurate angle calculations in collision cascades
- **Is the Second Program**: Indicates this is an improved version of a previous implementation
- **Laveraging Accurate Θ**: Emphasizes the program's focus on precise angle calculations in collision simulations

## Features

- Simulate collision cascades within materials.
- Compute irradiation damage (e.g., number of vacancies).
- Flexible geometry and customizable parameters.
- Supports a wide range of ion–material combinations.
- Provides convenient data output (CSV, dump files, etc.).

## Installation Requirements

- Julia 1.0 or higher.
- Dependencies:
  - LinearAlgebra
  - PeriodicTable
  - QuadGK
  - ProgressMeter (for displaying progress bars, configured in the code).

## Getting Started

1. Clone the repository to your local machine.
2. Ensure you have the required Julia packages installed.
3. Modify your calculation parameters and material settings as needed.
4. Run your simulation script.

## Usage

1. Define your material structure (lattice vectors, bounding box, basis atoms, etc.).
2. Assign simulation parameters (ion types, energies, collision parameters, etc.).
3. Initialize the simulator.
4. Execute your irradiation script.
5. Analyze the results (e.g., vacancy counts).

## Example

Below is a simplified example of simulating xenon irradiation on graphene:

```julia
# Load the BCA.jl modules
include("path/to/bca.jl")

# Define lattice parameters
a = 1.42
b = 3.35
primaryVectors = [3.0*a 0.0 0.0; 0.0 sqrt(3)*a 0.0; 0.0 0.0 b]
boxSizes = [10, 20, 10]
inputGridVectors = [2.1*a 0.0 0.0; 0.0 2.1*a 0.0; 0.0 0.0 2.1*a]
periodic = [true, true, false]
latticeRanges = [0 10; 0 20; 5 6]   
basis = [0.0 0.0 0.0; 1/3 0.0 0.0; 1/2 1/2 0.0; 5/6 1/2 0.0]
basisTypes = [1, 1, 1, 1]

# Simulation parameters
θτFileName = "CNe.thetatau"
pMax = 1.45
vacancyRecoverDistance = 1.3
typeDict = Dict(
    1 => Element("C", 22.0, 22.0),  
    2 => Element("Xe", 22.0, 22.0)  
)
# Default values in Parameters (for reference):
#   stopEnergy = 0.1
#   pLMax = 2.0
#   isDumpInCascade = false
#   isLog = false

parameters = Parameters(θτFileName, pMax, vacancyRecoverDistance, typeDict)

# Initialize the simulator
simulator = Simulator(primaryVectors, boxSizes, inputGridVectors, periodic,
    latticeRanges, basis, basisTypes, parameters)
Save!(simulator)

function Irradiation(simulator, energy)
    Restore!(simulator)
    simulator.nIrradiation += 1
    ionPosition = RandomInAnUnitGrapheneCell(1.39667) .+ [18.855045, 22.981482313368623, 20]
    ion = Atom(2, ionPosition, parameters)
    SetVelocityDirection!(ion, [0.0, 0.0, -1.0])
    SetEnergy!(ion, energy)
    push!(simulator, ion)
    Cascade!(ion, simulator)
    return CountVacancies(simulator)
end

computerNumberPerEnergy = 1000
vacancy_data = Dict{Int64, Vector{Int64}}()

energies = [10*1.2^x for x in 0:35]
for (n, energy) in enumerate(energies)
    vacancyNumbers = Vector{Int64}(undef, computerNumberPerEnergy)
    for i in 1:computerNumberPerEnergy
        vacancyNumbers[i] = Irradiation(simulator, energy)
    end
    vacancy_data[n] = vacancyNumbers
end

# Save results
open("vacancy.csv", "w") do f
    write(f, "n,energy,nVacancies\n")
    for n in sort(collect(keys(vacancy_data)))
        for vNum in vacancy_data[n]
            write(f, "$(n),$(energies[n]),$(vNum)\n")
        end
    end
end
```
  
## Parameters

- **Lattice / Box**:
  - `primaryVectors`: defines the unit cell of the material.
  - `boxSizes`: number of cells in each dimension.
  - `periodic`: boundary conditions (true or false).

- **Collision Settings**:
  - `θτFileName`: file with angle and displacement data.
  - `pMax`: maximum impact parameter.
  - `vacancyRecoverDistance`: threshold distance for vacancy recovery.
  - `stopEnergy` (default = 0.1 eV): minimum energy below which an atom stops moving.
  - `pLMax` (default = 2.0): maximum collision distance factor.
  - `isDumpInCascade` (default = false): whether to dump data after each collision event.
  - `isLog` (default = false): whether to print log info.

- **Ion/Material**:
  - `typeDict`: dictionary defining the element (name, mass, radius, etc.).

## Output

- `vacancy.csv`: CSV logging the energy and the corresponding vacancy count. 
- Optional dump files (`Cascade_*.dump`): track the intermediate steps of collision cascades.

Feel free to tweak the code as needed for your specific research or simulation goals. Happy simulating! ☀️

## References

[1] I2DM: A Monte Carlo framework for ion irradiation on two-dimensional materials
