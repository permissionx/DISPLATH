using LinearAlgebra
using StaticArrays
import Base: push!
using StableRNGs, Random, Base.Threads
using QuadGK
import Base: push!
import Base: delete!
using Interpolations   
using Dates
using ProgressMeter
using Distributions


#using PyCall
# @pyimport dscribe.descriptors as descriptors
#include("debug.jl")
#using .Recorder
include("debug.jl")
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

