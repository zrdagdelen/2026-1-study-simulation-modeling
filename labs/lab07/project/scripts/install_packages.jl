using DrWatson
@quickactivate "project"

using Pkg
Pkg.add.([
    "ConcurrentSim", "ResumableFunctions", "StableRNGs"
])
